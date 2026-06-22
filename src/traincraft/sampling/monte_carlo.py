"""Sampler: Metropolis Monte Carlo.

Moves:
  translate  — rigid translation of one fragment
  rotate     — rigid rotation about fragment centroid
  conformer  — RDKit ETKDG conformer swap, aligned to current placement

Fragment identity comes from tc_fragment array (set by builders/sources).
Falls back to treating the whole system as one rigid fragment if unset.
"""

from __future__ import annotations

import logging

import numpy as np
from ase import units
from ase.io.trajectory import Trajectory

from ..core import (
    Job,
    Provenance,
    Structure,
    fragment_ids,
    get_fragments,
    register,
    set_fragments,
)
from ..core.fragments import FRAMEWORK, fragment_mask, infer_fragments  # noqa: F401

logger = logging.getLogger(__name__)


def _fragment_centroid(atoms, mask: np.ndarray) -> np.ndarray:
    return atoms.get_positions()[mask].mean(axis=0)


def _move_translate(atoms, fid: int, rng: np.random.Generator, max_d: float) -> None:
    mask = fragment_mask(atoms, fid)
    delta = rng.uniform(-max_d, max_d, 3)
    pos = atoms.get_positions()
    pos[mask] += delta
    atoms.set_positions(pos)


def _move_rotate(atoms, fid: int, rng: np.random.Generator, max_deg: float) -> None:
    from scipy.spatial.transform import Rotation

    mask = fragment_mask(atoms, fid)
    pos = atoms.get_positions()
    centroid = pos[mask].mean(axis=0)

    axis = rng.standard_normal(3)
    axis /= np.linalg.norm(axis)
    angle_deg = rng.uniform(-max_deg, max_deg)
    rot = Rotation.from_rotvec(np.deg2rad(angle_deg) * axis)

    pos[mask] = rot.apply(pos[mask] - centroid) + centroid
    atoms.set_positions(pos)


def _move_conformer(
    atoms,
    fid: int,
    rng: np.random.Generator,
    fragment_smiles: dict[str, str] | None,
) -> bool:
    """Replace fragment atoms with a fresh ETKDG conformer aligned by centroid + PCA axes.

    Returns True if the move was applied, False if it was skipped.
    """
    if fragment_smiles is None:
        return False
    smiles = fragment_smiles.get(str(fid))
    if smiles is None:
        return False

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return False

    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        return False
    mol = Chem.AddHs(parsed)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(rng.integers(0, 2**31))
    if AllChem.EmbedMolecule(mol, params) < 0:
        return False
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer(0)
    new_pos = np.array(conf.GetPositions())
    new_syms = [a.GetSymbol() for a in mol.GetAtoms()]

    mask = fragment_mask(atoms, fid)
    old_pos = atoms.get_positions()[mask]

    # The new conformer's atoms are copied into the fragment slots index-for-index
    # below, so they MUST line up 1:1 — same count, same elements, same order — as
    # the fragment's atoms. If they don't (e.g. a canonical-SMILES rebuild reorders
    # the atoms), copying would drop each element's coordinates onto the WRONG atom
    # and silently corrupt the molecule (an OH's H landing on a CH3, etc.). Refuse
    # the move rather than emit a broken geometry.
    frag_syms = [s for s, m in zip(atoms.get_chemical_symbols(), mask, strict=True) if m]
    if new_syms != frag_syms:
        logger.warning(
            "monte_carlo: conformer move skipped — rebuilt atom order %s does not "
            "match fragment %d (%s); not applying to avoid corrupting the molecule",
            "".join(new_syms), fid, "".join(frag_syms),
        )
        return False

    old_centroid = old_pos.mean(axis=0)
    new_centroid = new_pos.mean(axis=0)
    new_pos = new_pos - new_centroid + old_centroid

    # Align principal axes via SVD (Kabsch-style).
    # Input shape (n_atoms, 3); full_matrices=False gives vh shape (3, 3).
    _, _, vh_old = np.linalg.svd(old_pos - old_centroid, full_matrices=False)
    _, _, vh_new = np.linalg.svd(new_pos - old_centroid, full_matrices=False)
    rot = vh_old.T @ vh_new
    # Ensure proper rotation (det=+1), not a reflection.
    if np.linalg.det(rot) < 0:
        vh_new[-1] *= -1
        rot = vh_old.T @ vh_new
    new_pos = (rot @ (new_pos - old_centroid).T).T + old_centroid

    pos = atoms.get_positions()
    pos[mask] = new_pos
    atoms.set_positions(pos)
    return True


@register("sampler", "monte_carlo")
def sample_monte_carlo(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    atoms = structure.atoms.copy()
    atoms.calc = calc

    rng = np.random.default_rng(cfg.seed)

    # Ensure fragments are set.
    if get_fragments(atoms) is None:
        logger.warning(
            "monte_carlo: no tc_fragment array found; "
            "treating entire system as one rigid fragment"
        )
        set_fragments(atoms, np.zeros(len(atoms), dtype=int))

    mobile = fragment_ids(atoms)
    if not mobile:
        raise ValueError(
            "monte_carlo: no mobile fragments to move "
            "(all atoms are framework or no fragment array set)"
        )

    # Normalise move probabilities.
    raw = np.array([cfg.p_translate, cfg.p_rotate, cfg.p_conformer], dtype=float)
    probs = raw / raw.sum()

    fragment_smiles: dict[str, str] | None = structure.provenance.extra.get(
        "fragment_smiles"
    )
    conformer_warned = False

    e_current = atoms.get_potential_energy()
    kbt = units.kB * cfg.temperature

    frames: list[Structure] = []
    traj = Trajectory(str(job.path("mc.traj")), "w", atoms)

    n_accepted = 0

    for step in range(cfg.steps):
        fid = int(rng.choice(mobile))

        trial = atoms.copy()
        trial.calc = calc

        move = rng.choice(["translate", "rotate", "conformer"], p=probs)

        if move == "translate":
            _move_translate(trial, fid, rng, cfg.max_translate)
        elif move == "rotate":
            _move_rotate(trial, fid, rng, cfg.max_rotate)
        else:
            applied = _move_conformer(trial, fid, rng, fragment_smiles)
            if not applied:
                if not conformer_warned:
                    logger.warning(
                        "monte_carlo: conformer move skipped "
                        "(rdkit unavailable or no SMILES recorded for fragment %d); "
                        "falling back to translate",
                        fid,
                    )
                    conformer_warned = True
                _move_translate(trial, fid, rng, cfg.max_translate)

        e_trial = trial.get_potential_energy()
        delta_e = e_trial - e_current

        accepted = bool(delta_e <= 0 or rng.random() < np.exp(-delta_e / kbt))
        if accepted:
            atoms = trial
            e_current = e_trial
            n_accepted += 1

            if cfg.refresh_fragments:
                new_frag = infer_fragments(atoms, scale=cfg.refresh_scale)
                old_count = len(mobile)
                set_fragments(atoms, new_frag)
                mobile = fragment_ids(atoms)
                if len(mobile) != old_count:
                    logger.info(
                        "monte_carlo step %d: fragment count changed %d → %d",
                        step, old_count, len(mobile),
                    )

        if (step + 1) % cfg.interval == 0:
            # "accepted": snapshot the chain state (Boltzmann ensemble). "trials":
            # snapshot the *proposed* geometry whether accepted or not, so the
            # high-energy, off-equilibrium configs a Metropolis filter would discard
            # are kept — broad coverage for an unbiased MLIP training set. The trial's
            # energy is already computed above; only its forces are evaluated here.
            snapshot = trial if cfg.record == "trials" else atoms
            props: dict = {}
            try:
                props["energy"] = float(snapshot.get_potential_energy())
                props["forces"] = snapshot.get_forces().copy()
            except Exception:  # noqa: BLE001
                pass

            frames.append(
                Structure.from_ase(
                    snapshot,
                    properties=props,
                    provenance=Provenance(
                        origin="ml_sampled",
                        source=f"sampler:monte_carlo:{cfg.record}",
                        calculator=getattr(calc, "name", calc.__class__.__name__),
                        parents=[structure.hash],
                        extra={"mc_accepted": accepted},
                    ),
                )
            )
            traj.write(snapshot)

    logger.info(
        "MC finished: %d steps, %d accepted (%.1f%%), %d frames recorded",
        cfg.steps, n_accepted, 100 * n_accepted / max(cfg.steps, 1), len(frames),
    )
    return frames
