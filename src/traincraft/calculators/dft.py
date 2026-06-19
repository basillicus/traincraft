"""DFT labeling calculators: FHI-aims and Quantum ESPRESSO.

These are the *labelers* of the pipeline (DESIGN §6, §8): they produce
energy/forces/stress for every frame, and on request dipole (SCF) and
polarizability (DFPT — linear response, materially heavier than SCF).

The plugins are deliberately **container-agnostic** (DESIGN §20.3): the run
command is never hard-coded. It is read from an environment variable so that an
HPC executor can inject ``srun apptainer exec … traincraft-dft.sif aims.x``
while local/dev just runs a bare binary. ``core`` writes the inputs and parses
the outputs; the ``.sif`` is a pure DFT worker.

Heavy imports (ASE FileIO calculators) are deferred inside the factory bodies so
the package imports with no DFT stack present, mirroring ``potentials.py``.
"""

from __future__ import annotations

import os
from pathlib import PurePath

from ..core import register

# energy/forces/stress are produced by every SCF; dipole/polarizability are opt-in.
_EFS = {"energy", "forces", "stress"}
_DFT_CAPS = _EFS | {"dipole", "polarizability"}

# Environment variables for command injection (see DESIGN §20.3).
AIMS_COMMAND_ENV = "TRAINCRAFT_AIMS_COMMAND"
AIMS_SPECIES_ENV = ("TRAINCRAFT_AIMS_SPECIES_DIR", "AIMS_SPECIES_DIR")
PW_COMMAND_ENV = "TRAINCRAFT_PW_COMMAND"
PW_PSEUDO_ENV = ("TRAINCRAFT_PW_PSEUDO_DIR", "ESPRESSO_PSEUDO")

DEFAULT_AIMS_COMMAND = "aims.x"
DEFAULT_PW_COMMAND = "pw.x"


def _env_first(names: tuple[str, ...]) -> str | None:
    """Return the first set, non-empty environment variable from ``names``."""
    for name in names:
        value = os.environ.get(name)
        if value:
            return value
    return None


@register("calculator", "fhi_aims", capabilities=_DFT_CAPS)
def build_fhi_aims(cfg):
    """ASE FHI-aims ``GenericFileIOCalculator``.

    Capabilities: energy/forces/stress always; dipole and polarizability when
    listed in ``cfg.properties``. The run command comes from
    ``$TRAINCRAFT_AIMS_COMMAND`` (default ``aims.x``); the species directory from
    ``cfg.species_dir`` or ``$TRAINCRAFT_AIMS_SPECIES_DIR`` / ``$AIMS_SPECIES_DIR``.

    Polarizability is requested via DFPT. The control.in keyword differs by
    boundary conditions: ``DFPT dielectric`` for periodic systems, ``DFPT
    polarizability`` for molecules. Whether the system is periodic is taken from
    ``cfg.periodic`` (auto-True when a k-grid is supplied).
    """
    from ase.calculators.aims import Aims, AimsProfile

    command = os.environ.get(AIMS_COMMAND_ENV) or DEFAULT_AIMS_COMMAND

    # The basis-set level (light/tight/...) is selected by the species *path*,
    # not a control.in keyword. If the configured dir is a defaults root, append
    # the level; otherwise assume it already points at the level directory.
    species_root = cfg.species_dir or _env_first(AIMS_SPECIES_ENV)
    species_dir = species_root
    if species_root is not None and cfg.species_defaults:
        if PurePath(species_root).name != cfg.species_defaults:
            species_dir = str(PurePath(species_root) / cfg.species_defaults)

    profile = AimsProfile(command=command, default_species_directory=species_dir)

    requested = set(cfg.properties)
    periodic = cfg.periodic if cfg.periodic is not None else cfg.kpts is not None

    # ASE-Aims standard + native control.in keywords. (species_defaults is NOT a
    # control.in keyword — it is folded into species_dir above.)
    parameters: dict = {
        "xc": cfg.xc,
        "species_dir": species_dir,
        "relativistic": cfg.relativistic,
    }
    if species_dir is None:
        # Let the profile's default_species_directory drive it instead.
        parameters.pop("species_dir")
    if cfg.spin != "none":
        parameters["spin"] = cfg.spin
    if cfg.kpts is not None:
        parameters["k_grid"] = tuple(cfg.kpts)

    # dipole: control.in `output dipole`. ASE-Aims appends 'dipole' to `output`
    # automatically when 'dipole' is in the requested properties at calc time,
    # but we also set it explicitly so a bare label run emits it.
    if "dipole" in requested:
        parameters.setdefault("output", [])
        outputs = list(parameters["output"])
        if "dipole" not in outputs:
            outputs.append("dipole")
        parameters["output"] = outputs

    # polarizability via DFPT (linear response).
    if "polarizability" in requested:
        parameters["dfpt"] = "dielectric" if periodic else "polarizability"

    # passthrough: arbitrary control.in keywords win over our defaults.
    parameters.update(cfg.extra)

    calc = Aims(profile=profile, **parameters)
    return calc


@register("calculator", "qe", capabilities=_EFS | {"dipole"})
def build_qe(cfg):
    """ASE Quantum ESPRESSO (``pw.x``) calculator.

    Capabilities: energy/forces/stress always; dipole via SCF (Berry-phase /
    ``dipfield``). The command comes from ``$TRAINCRAFT_PW_COMMAND`` (default
    ``pw.x``); pseudopotentials directory from ``cfg.pseudo_dir`` or
    ``$TRAINCRAFT_PW_PSEUDO_DIR`` / ``$ESPRESSO_PSEUDO``.

    Polarizability is **not** wired here: in QE it requires a separate ``ph.x``
    (DFPT) run after the SCF, which is a multi-binary workflow outside the scope
    of a single ASE ``FileIOCalculator``. Requesting it raises with a clear
    message; the extension point is :func:`build_qe` plus a future ``ph.x``
    profile. Hence ``polarizability`` is intentionally absent from the declared
    capabilities for this calculator.
    """
    from ase.calculators.espresso import Espresso, EspressoProfile

    requested = set(cfg.properties)
    if "polarizability" in requested:
        raise NotImplementedError(
            "QE polarizability needs a separate ph.x (DFPT) run, which is not "
            "wired into the single-binary pw.x calculator. Use the 'fhi_aims' "
            "calculator for polarizability, or extend build_qe with a ph.x "
            "profile."
        )

    command = os.environ.get(PW_COMMAND_ENV) or DEFAULT_PW_COMMAND
    pseudo_dir = cfg.pseudo_dir or _env_first(PW_PSEUDO_ENV)
    if pseudo_dir is None:
        # EspressoProfile requires a pseudo_dir; keep an explicit, debuggable hint.
        pseudo_dir = "."

    profile = EspressoProfile(command=command, pseudo_dir=pseudo_dir)

    # Build the pw.x control/system namelists. input_data wins over our defaults.
    input_data: dict = {
        "control": {"calculation": "scf", "tprnfor": True, "tstress": True},
        "system": {"ecutwfc": cfg.ecutwfc},
    }
    if cfg.ecutrho is not None:
        input_data["system"]["ecutrho"] = cfg.ecutrho
    if "dipole" in requested:
        # Berry-phase dipole correction; user can refine via cfg.input_data.
        input_data["control"]["dipfield"] = True
        input_data["control"].setdefault("tefield", False)

    # merge user-supplied namelists (nested dict) on top of defaults
    for section, values in cfg.input_data.items():
        if isinstance(values, dict):
            input_data.setdefault(section, {}).update(values)
        else:
            input_data[section] = values

    kwargs: dict = {
        "pseudopotentials": dict(cfg.pseudopotentials),
        "input_data": input_data,
    }
    if cfg.kspacing is not None:
        kwargs["kspacing"] = cfg.kspacing
    elif cfg.kpts is not None:
        kwargs["kpts"] = tuple(cfg.kpts)

    calc = Espresso(profile=profile, **kwargs)
    return calc
