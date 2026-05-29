"""Sampler: Langevin MD (ported from legacy ``samplers.sampling_from_MD``).

Parameterized (no global config), writes into the job dir, returns frames with
provenance. No ``os.chdir``, no uuid filename side-channel.
"""

from __future__ import annotations

import logging

from ase import units
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from ..core import Job, Provenance, Structure, register

logger = logging.getLogger(__name__)


@register("sampler", "md")
def sample_md(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    atoms = structure.atoms.copy()
    atoms.calc = calc
    MaxwellBoltzmannDistribution(atoms, temperature_K=cfg.temperature)

    dyn = Langevin(
        atoms,
        timestep=cfg.timestep * units.fs,
        temperature_K=cfg.temperature,
        friction=cfg.friction,
    )

    frames: list[Structure] = []

    def _grab() -> None:
        props: dict = {}
        try:
            props["energy"] = float(atoms.get_potential_energy())
            props["forces"] = atoms.get_forces().copy()
        except Exception:  # noqa: BLE001 - properties are best-effort here
            pass
        frames.append(
            Structure.from_ase(
                atoms,
                properties=props,
                provenance=Provenance(
                    origin="ml_sampled",
                    source="sampler:md",
                    calculator=getattr(calc, "name", calc.__class__.__name__),
                    parents=[structure.hash],
                ),
            )
        )

    traj = Trajectory(str(job.path("md.traj")), "w", atoms)
    dyn.attach(traj.write, interval=cfg.interval)
    dyn.attach(_grab, interval=cfg.interval)

    logger.info("MD: T=%s steps=%s interval=%s", cfg.temperature, cfg.steps, cfg.interval)
    dyn.run(cfg.steps)
    logger.info("MD finished: %d frames sampled", len(frames))
    return frames
