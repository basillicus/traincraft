"""Tests for the DFT labeling calculators (fhi_aims, qe).

The heavy DFT binaries are not installed, so nothing is *run*: we only build the
ASE calculators and inspect their parameters/profile. ASE is a core dep.
"""

from __future__ import annotations

import pytest

from traincraft.calculators import make_calculator
from traincraft.calculators.dft import (
    AIMS_COMMAND_ENV,
    AIMS_SPECIES_ENV,
    DEFAULT_AIMS_COMMAND,
    DEFAULT_PW_COMMAND,
    PW_COMMAND_ENV,
)
from traincraft.config.models import FhiAimsCalc, QeCalc
from traincraft.core import RegistryError, capabilities, get

pytest.importorskip("ase")


# --------------------------------------------------------------------------- registry
def test_factories_registered():
    assert get("calculator", "fhi_aims") is not None
    assert get("calculator", "qe") is not None


def test_capabilities():
    aims_caps = capabilities("calculator", "fhi_aims")
    assert {"energy", "forces", "stress", "dipole", "polarizability"} <= aims_caps
    qe_caps = capabilities("calculator", "qe")
    assert {"energy", "forces", "stress", "dipole"} <= qe_caps
    # QE does not advertise polarizability (needs ph.x)
    assert "polarizability" not in qe_caps


# ----------------------------------------------------------------------- aims build
def test_aims_resolvable_via_make_calculator():
    calc = make_calculator(FhiAimsCalc(xc="pbe", species_dir="/tmp/species/light"))
    assert type(calc).__name__ == "Aims"
    assert calc.parameters["xc"] == "pbe"


def test_aims_default_command(monkeypatch):
    monkeypatch.delenv(AIMS_COMMAND_ENV, raising=False)
    calc = make_calculator(FhiAimsCalc(species_dir="/tmp/light"))
    assert calc.profile.command == DEFAULT_AIMS_COMMAND


def test_aims_command_injection(monkeypatch):
    injected = "srun apptainer exec --bind /scratch traincraft-dft.sif aims.x"
    monkeypatch.setenv(AIMS_COMMAND_ENV, injected)
    calc = make_calculator(FhiAimsCalc(species_dir="/tmp/light"))
    assert calc.profile.command == injected
    # the profile splits it into argv when running
    assert calc.profile._split_command[0] == "srun"


def test_aims_species_dir_from_env(monkeypatch):
    for var in AIMS_SPECIES_ENV:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv(AIMS_SPECIES_ENV[0], "/opt/species/defaults_2020")
    calc = make_calculator(FhiAimsCalc(species_defaults="tight"))
    # species_defaults folded onto the env-provided root
    assert calc.profile.default_species_directory == "/opt/species/defaults_2020/tight"
    assert calc.parameters["species_dir"] == "/opt/species/defaults_2020/tight"


def test_aims_species_dir_already_level():
    # if the dir already points at the level, don't double-append
    calc = make_calculator(FhiAimsCalc(species_dir="/opt/defaults_2020/light",
                                       species_defaults="light"))
    assert calc.profile.default_species_directory == "/opt/defaults_2020/light"


def test_aims_dipole_output():
    calc = make_calculator(
        FhiAimsCalc(species_dir="/tmp/light", properties=["dipole"])
    )
    assert "dipole" in calc.parameters["output"]


def test_aims_polarizability_molecular_uses_dfpt_polarizability():
    # non-periodic (no kpts, periodic unset) -> molecular DFPT polarizability
    calc = make_calculator(
        FhiAimsCalc(species_dir="/tmp/light", properties=["polarizability"])
    )
    assert calc.parameters["dfpt"] == "polarizability"


def test_aims_polarizability_periodic_uses_dfpt_dielectric():
    # k-grid implies periodic -> DFPT dielectric
    calc = make_calculator(
        FhiAimsCalc(
            species_dir="/tmp/light",
            properties=["polarizability"],
            kpts=(4, 4, 4),
        )
    )
    assert calc.parameters["dfpt"] == "dielectric"
    assert calc.parameters["k_grid"] == (4, 4, 4)


def test_aims_explicit_periodic_flag_overrides():
    # explicit periodic=True without kpts still selects dielectric
    calc = make_calculator(
        FhiAimsCalc(
            species_dir="/tmp/light",
            properties=["polarizability"],
            periodic=True,
        )
    )
    assert calc.parameters["dfpt"] == "dielectric"


def test_aims_no_polarizability_no_dfpt():
    calc = make_calculator(FhiAimsCalc(species_dir="/tmp/light"))
    assert "dfpt" not in calc.parameters


def test_aims_extra_passthrough_wins():
    calc = make_calculator(
        FhiAimsCalc(species_dir="/tmp/light", xc="pbe", extra={"xc": "pbe0",
                                                               "sc_accuracy_rho": 1e-5})
    )
    assert calc.parameters["xc"] == "pbe0"
    assert calc.parameters["sc_accuracy_rho"] == 1e-5


def test_aims_species_defaults_not_a_control_keyword():
    # species_defaults must never leak into control.in parameters
    calc = make_calculator(FhiAimsCalc(species_dir="/tmp/light",
                                       species_defaults="tight"))
    assert "species_defaults" not in calc.parameters


# ------------------------------------------------------------------------- qe build
def test_qe_resolvable_via_make_calculator():
    calc = make_calculator(QeCalc(pseudo_dir="/tmp/pseudo"))
    assert type(calc).__name__ == "Espresso"


def test_qe_default_command(monkeypatch):
    monkeypatch.delenv(PW_COMMAND_ENV, raising=False)
    calc = make_calculator(QeCalc(pseudo_dir="/tmp/pseudo"))
    assert calc.profile.command == DEFAULT_PW_COMMAND


def test_qe_command_injection(monkeypatch):
    injected = "srun apptainer exec traincraft-dft.sif pw.x"
    monkeypatch.setenv(PW_COMMAND_ENV, injected)
    calc = make_calculator(QeCalc(pseudo_dir="/tmp/pseudo"))
    assert calc.profile.command == injected


def test_qe_pseudo_dir_and_inputs():
    calc = make_calculator(
        QeCalc(
            pseudo_dir="/tmp/pseudo",
            ecutwfc=80.0,
            kpts=(2, 2, 2),
            pseudopotentials={"Si": "Si.upf"},
        )
    )
    assert calc.profile.pseudo_dir == "/tmp/pseudo"
    assert calc.parameters["input_data"]["system"]["ecutwfc"] == 80.0
    assert calc.parameters["input_data"]["control"]["calculation"] == "scf"
    assert calc.parameters["kpts"] == (2, 2, 2)
    assert calc.parameters["pseudopotentials"] == {"Si": "Si.upf"}


def test_qe_dipole_sets_dipfield():
    calc = make_calculator(QeCalc(pseudo_dir="/tmp/pseudo", properties=["dipole"]))
    assert calc.parameters["input_data"]["control"]["dipfield"] is True


def test_qe_input_data_override():
    calc = make_calculator(
        QeCalc(
            pseudo_dir="/tmp/pseudo",
            input_data={"system": {"occupations": "smearing", "degauss": 0.01}},
        )
    )
    system = calc.parameters["input_data"]["system"]
    assert system["occupations"] == "smearing"
    assert system["degauss"] == 0.01
    # defaults preserved
    assert system["ecutwfc"] == 60.0


def test_qe_polarizability_raises():
    with pytest.raises(NotImplementedError, match="ph.x"):
        make_calculator(QeCalc(pseudo_dir="/tmp/pseudo",
                               properties=["polarizability"]))


# ----------------------------------------------------------------- config guardrails
def test_unknown_dft_calculator_name():
    with pytest.raises(RegistryError):
        get("calculator", "vasp")


def test_dft_config_forbids_extra_fields():
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        FhiAimsCalc(species_dir="/tmp", bogus_field=1)
