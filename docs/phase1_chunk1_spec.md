# Phase 1 — Chunk 1 Spec: Molecules-on-Surfaces + Monte Carlo

**Status:** ready to implement
**Audience:** an implementer (small model or human) who follows this literally.
**Rule:** no architecture decisions are left open here. If something is ambiguous,
it is a bug in this spec — stop and ask, do not improvise.

This chunk delivers the hard abstraction of Phase 1 — *fragment identity* — and the
geometry + sampling features that depend on it: a SMILES source, two molecule-on-
surface builders (single adsorbate via ASE, multi-molecule coverage via Packmol),
and the Metropolis Monte Carlo sampler (rigid-body moves + RDKit conformer swaps).
The mechanical crystal/defect/slab builders are a **later** chunk and are out of scope.

---

## 0. Design decisions (settled — do not relitigate)

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | **Fragment identity is a per-atom integer array** `tc_fragment` in `atoms.arrays`: `-1` = framework/substrate (immobile), `0,1,2,…` = mobile fragment ids. | Survives `Atoms.copy()`, extxyz round-trip, supercell `repeat()`, and MD trajectory writing **for free**, because ASE propagates `arrays`. No `Structure` schema change. |
| D2 | **Runtime connectivity inference is also provided**, as an explicit helper `infer_fragments(atoms, scale=1.2)`, NOT as the default. It is used (a) by builders that need to tag an assembled system, and (b) opt-in by the MC sampler via `refresh_fragments` for reactive runs where bonds change. | The user wants both: the stored array is simple and right for initial geometry manipulation; inference is the escape hatch for chemistry that changes connectivity mid-run. |
| D3 | **Both placement mechanisms ship now:** `surface_adsorbate` builder (single molecule, `ase.build.add_adsorbate`, zero extra deps) and `surface_packing` builder (N molecules, Packmol via conda-forge). | User explicitly chose "both from the start." |
| D4 | RDKit and Packmol are **optional deps**, imported lazily inside the function that needs them, with a clear `ImportError` message naming the pixi env. Mirrors the existing `rattle.py`/hiphive pattern. | Keeps the `default` env light; `science` env carries them. |
| D5 | The MC sampler is **molecule-aware via `tc_fragment`**. If the array is absent it falls back to treating the whole non-framework system as a single rigid fragment, and logs a warning. It never silently does nothing. | A sampler that needs fragments must degrade predictably. |
| D6 | Fragment tagging is a **first-class concern of `Structure`**: add thin helpers `Structure.fragments` (read), `Structure.set_fragments(arr)` (write), `Structure.n_fragments`. These wrap `atoms.arrays["tc_fragment"]`; they do not add dataclass fields. | One obvious API; keeps the per-atom-array decision (D1) encapsulated so callers never touch the raw key. |

---

## 1. The fragment-identity layer (build FIRST — everything depends on it)

### 1.1 `src/traincraft/core/fragments.py` (new)

The single source of truth for what a "fragment" means. Pure functions, no ASE
calculator, no IO.

```python
"""Fragment identity: which atoms move together.

Convention (project-wide):
  tc_fragment[i] == -1   -> atom i is framework / substrate (immobile)
  tc_fragment[i] >= 0    -> atom i belongs to mobile fragment with that id

The array lives in atoms.arrays["tc_fragment"] so it survives copy / extxyz /
MD / supercell automatically (ASE propagates per-atom arrays).
"""

from __future__ import annotations

import numpy as np
from ase import Atoms

FRAGMENT_KEY = "tc_fragment"
FRAMEWORK = -1


def get_fragments(atoms: Atoms) -> np.ndarray | None:
    """Return the per-atom fragment array, or None if unset."""
    if FRAGMENT_KEY in atoms.arrays:
        return atoms.arrays[FRAGMENT_KEY].astype(int)
    return None


def set_fragments(atoms: Atoms, frag: np.ndarray | list[int]) -> None:
    """Attach/overwrite the per-atom fragment array (length must equal len(atoms))."""
    frag = np.asarray(frag, dtype=int)
    if frag.shape != (len(atoms),):
        raise ValueError(f"fragment array must have shape ({len(atoms)},), got {frag.shape}")
    atoms.set_array(FRAGMENT_KEY, frag)


def fragment_ids(atoms: Atoms) -> list[int]:
    """Sorted list of mobile fragment ids (excludes FRAMEWORK)."""
    frag = get_fragments(atoms)
    if frag is None:
        return []
    return sorted(int(i) for i in np.unique(frag) if i != FRAMEWORK)


def fragment_mask(atoms: Atoms, fid: int) -> np.ndarray:
    """Boolean mask selecting atoms of fragment `fid`."""
    frag = get_fragments(atoms)
    if frag is None:
        raise ValueError("no fragment array set on these atoms")
    return frag == fid


def infer_fragments(atoms: Atoms, scale: float = 1.2, framework_mask=None) -> np.ndarray:
    """Assign fragment ids by connected components of a covalent-radius graph.

    Two atoms bond if distance < scale * (r_cov[i] + r_cov[j]), using ASE
    natural_cutoffs + NeighborList with self_interaction=False, bothways=True.
    `framework_mask` (optional bool array): atoms marked True are forced to
    FRAMEWORK and excluded from the connectivity graph. Returns the array; does
    NOT mutate `atoms` (caller decides whether to set_fragments).

    Implementation: ase.neighborlist.natural_cutoffs + NeighborList +
    scipy.sparse.csgraph.connected_components on the adjacency. scipy is already
    an ASE dependency.
    """
    ...
```

**Acceptance for this file:**
- `set_fragments` rejects wrong-length arrays.
- `infer_fragments` on an isolated H2O + CO2 (placed far apart, no cell) returns 2
  components → ids `{0, 1}`.
- `infer_fragments` with a `framework_mask` over a slab+molecule returns the slab
  atoms as `-1` and the molecule as one component.
- Round-trip: `set_fragments` → `write_frames` → `read_frames` →
  `get_fragments` returns the same array. **(This requires §1.3.)**

### 1.2 `Structure` helpers — `src/traincraft/core/structure.py` (edit)

Add three thin properties/methods that delegate to `fragments.py`. Do **not** add
dataclass fields.

```python
# in Structure:
@property
def fragments(self):
    from .fragments import get_fragments
    return get_fragments(self.atoms)

def set_fragments(self, frag) -> None:
    from .fragments import set_fragments
    set_fragments(self.atoms, frag)

@property
def n_fragments(self) -> int:
    from .fragments import fragment_ids
    return len(fragment_ids(self.atoms))
```

Export `get_fragments, set_fragments, infer_fragments, FRAGMENT_KEY, FRAMEWORK` from
`core/__init__.py` alongside the existing exports.

### 1.3 Make `tc_fragment` survive extxyz IO — `datasets/io.py` (edit)

`tc_fragment` lives in `atoms.arrays`, but `read_frames` only currently lifts
`tc_forces` out of arrays back into `properties`. Fragment is **not** a property —
it must stay in `atoms.arrays` on the reconstructed `Structure.atoms`. ASE's extxyz
writer/reader round-trips named per-atom arrays automatically, so the only required
change is a **test** confirming it (no code change expected). Add an assertion in the
io test (§7) that `tc_fragment` survives. If, and only if, ASE drops the column,
add it to the explicit columns list in `write` via `columns=[...]`. Verify first.

---

## 2. Config models — `src/traincraft/config/models.py` (edit)

Add the new source, two builders, and extend the MC sampling model. All extend
`TCModel` (so `extra="forbid"` is inherited). Add each to its discriminated union.

### 2.1 New source: SMILES

```python
class SmilesSource(TCModel):
    type: Literal["smiles"] = "smiles"
    smiles: str                         # e.g. "O=C=O" for CO2
    n_conformers: int = 1               # ETKDG conformers to generate; first is returned
    optimize: bool = True               # MMFF94 optimize after embedding
    seed: int | None = None             # ETKDG random seed (falls back to run.seed)
    vacuum: float = 6.0                 # box padding for the isolated molecule
```
Add to: `SourceConfig = Annotated[FileSource | ScratchSource | SmilesSource, Field(discriminator="type")]`

### 2.2 New builder: single adsorbate on a surface (ASE, zero extra deps)

```python
class SurfaceAdsorbateBuilder(TCModel):
    type: Literal["surface_adsorbate"] = "surface_adsorbate"
    # --- substrate (built with ase.build.fcc111/bcc110/hcp0001 family) ---
    element: str                        # e.g. "Cu"
    facet: Literal["fcc111","fcc100","fcc110","bcc110","bcc100","hcp0001"] = "fcc111"
    size: tuple[int, int, int] = (3, 3, 4)   # (nx, ny, layers)
    vacuum: float = 12.0                # total vacuum added along z
    # --- adsorbate: exactly one of molecule_name | smiles | file ---
    molecule_name: str | None = None    # ase.build.molecule g2 name (e.g. "CO")
    smiles: str | None = None           # RDKit-built adsorbate
    file: str | None = None             # path to an adsorbate geometry (resolved like FileSource)
    # --- placement ---
    site: Literal["ontop","bridge","hollow","fcc","hcp"] = "ontop"
    height: float = 2.0                 # Å above the surface
    offset: tuple[float, float] | None = None  # optional (x,y) fractional offset

    @model_validator(mode="after")
    def _one_adsorbate(self):
        n = sum(x is not None for x in (self.molecule_name, self.smiles, self.file))
        if n != 1:
            raise ValueError("surface_adsorbate needs exactly one of "
                             "'molecule_name', 'smiles', or 'file'")
        return self
```

### 2.3 New builder: multi-molecule coverage via Packmol

```python
class SurfacePackingBuilder(TCModel):
    type: Literal["surface_packing"] = "surface_packing"
    element: str
    facet: Literal["fcc111","fcc100","fcc110","bcc110","bcc100","hcp0001"] = "fcc111"
    size: tuple[int, int, int] = (4, 4, 4)
    vacuum: float = 20.0
    # adsorbate molecules to pack above the slab
    molecule_name: str | None = None
    smiles: str | None = None
    file: str | None = None
    n_molecules: int = 4
    tolerance: float = 2.0              # Packmol min interatomic distance (Å)
    region_height: float = 8.0          # thickness (Å) of the packing slab above the surface
    gap: float = 2.0                    # gap (Å) between surface top and packing region
    seed: int | None = None             # Packmol seed (falls back to run.seed)

    @model_validator(mode="after")
    def _one_adsorbate(self):
        n = sum(x is not None for x in (self.molecule_name, self.smiles, self.file))
        if n != 1:
            raise ValueError("surface_packing needs exactly one of "
                             "'molecule_name', 'smiles', or 'file'")
        return self
```

Add both to:
`BuilderConfig = Annotated[NanotubeBuilder | MoleculeBuilder | SurfaceAdsorbateBuilder | SurfacePackingBuilder, Field(discriminator="type")]`

### 2.4 Extend Monte Carlo sampling

Replace the existing `MonteCarloSampling` with:

```python
class MonteCarloSampling(TCModel):
    type: Literal["monte_carlo"] = "monte_carlo"
    steps: int = 1000
    temperature: float = 500.0          # Metropolis acceptance temperature (K)
    interval: int = 20                  # record a frame every N accepted-or-tried steps
    # move mix (probabilities; normalized internally; must not all be zero)
    p_translate: float = 0.5
    p_rotate: float = 0.4
    p_conformer: float = 0.1            # RDKit conformer swap (needs rdkit + fragment came from smiles)
    # move magnitudes
    max_translate: float = 0.5          # Å, max per-step rigid translation
    max_rotate: float = 30.0            # deg, max per-step rigid rotation
    # reactivity
    refresh_fragments: bool = False     # re-infer fragments each accepted step (D2)
    refresh_scale: float = 1.2          # passed to infer_fragments when refreshing
    seed: int | None = None

    @model_validator(mode="after")
    def _moves_nonzero(self):
        if self.p_translate + self.p_rotate + self.p_conformer <= 0:
            raise ValueError("monte_carlo needs at least one move probability > 0")
        return self
```

The existing `conformers: bool` field is **removed** (superseded by `p_conformer`).
Update example 06 / any config that referenced it — grep first.

---

## 3. SMILES source — `src/traincraft/geometry/sources/smiles.py` (new)

```python
@register("source", "smiles")
def source_smiles(cfg) -> Structure:
    # lazy import rdkit with helpful message naming the science env
    # 1. MolFromSmiles -> AddHs
    # 2. EmbedMultipleConfs (ETKDGv3), n=cfg.n_conformers, randomSeed=cfg.seed or 0xF00D
    # 3. if cfg.optimize: MMFFOptimizeMoleculeConfs
    # 4. take conformer 0 -> positions + atomic numbers -> ase.Atoms
    # 5. atoms.center(vacuum=cfg.vacuum/2); pbc stays False (isolated molecule)
    # 6. set_fragments(atoms, np.zeros(len(atoms)))   # whole molecule = fragment 0
    # provenance: origin="generated", source=f"source:smiles:{cfg.smiles}",
    #             extra={"smiles": cfg.smiles, "n_conformers": cfg.n_conformers}
    #   store the canonical smiles in provenance.extra so the MC conformer move
    #   can rebuild the molecule later (see §5.4).
```

**ImportError message (use verbatim style of rattle.py):**
> `"SMILES source needs RDKit. Install it with: pixi install -e science  (or: uv pip install rdkit)"`

Register it by importing in `geometry/sources/__init__.py` (match existing pattern).

**Acceptance:** `source_smiles` on `"O"` (water) returns a 3-atom `Structure`,
fragments all `0`, `provenance.extra["smiles"] == "O"`. Test skips if rdkit absent.

---

## 4. Surface builders — `src/traincraft/geometry/builders/surface.py` (new)

One file, both builders, plus a shared private adsorbate-resolver.

### 4.1 Shared helper

```python
def _resolve_adsorbate(molecule_name, smiles, file) -> Atoms:
    """Return an ase.Atoms for the adsorbate from whichever field is set.
    - molecule_name -> ase.build.molecule
    - smiles        -> reuse source_smiles logic (build a 1-conformer molecule)
    - file          -> ase.io.read (path already resolved by the loader; see note)
    """
```
Note: `file` paths inside builders are **not** auto-resolved by the existing
`_resolve_input_paths` (which only handles `[geometry.source].path`). Extend
`_resolve_input_paths` in `config/loader.py` to also resolve `builder.file` when
present (same relative-to-config rule). Add a test mirroring the existing
source-path test.

### 4.2 `surface_adsorbate`

```python
@register("builder", "surface_adsorbate")
def build_surface_adsorbate(cfg) -> Structure:
    # 1. slab = ase.build.{facet}(cfg.element, size=cfg.size, vacuum=cfg.vacuum/2)
    #    (map facet string -> the ase.build function: fcc111, fcc100, ... hcp0001)
    # 2. n_slab = len(slab)
    # 3. ads = _resolve_adsorbate(...)
    # 4. ase.build.add_adsorbate(slab, ads, height=cfg.height, position=cfg.site,
    #                            offset=cfg.offset)
    # 5. frag = full(-1) for the n_slab substrate atoms, 0 for the adsorbate atoms
    #    set_fragments(slab, frag)
    # provenance: origin="generated",
    #   source=f"builder:surface_adsorbate:{cfg.element}-{cfg.facet}"
    #   extra carries the adsorbate smiles if cfg.smiles set (for conformer moves)
```

### 4.3 `surface_packing`

```python
@register("builder", "surface_packing")
def build_surface_packing(cfg) -> Structure:
    # 1. build slab as above (no adsorbate)
    # 2. determine the packing box: x,y span the slab cell; z from
    #    (slab_top + gap) to (slab_top + gap + region_height)
    # 3. run Packmol (lazy import; see below) to place cfg.n_molecules copies of
    #    the adsorbate inside that box with cfg.tolerance
    # 4. concatenate slab + packed molecules
    # 5. fragments: slab atoms -> -1 ; each packed molecule -> 0,1,2,... by the
    #    known per-molecule atom count (n_atoms_adsorbate)
    # provenance as above, extra includes n_molecules and smiles if any
```

**Packmol invocation:** prefer calling the `packmol` binary (conda-forge) via a
temporary input deck written to `tempfile` (NOT the job dir — builders don't get a
Job), parsing the output `.xyz`. Lazy `shutil.which("packmol")` check;
if missing raise:
> `"surface_packing needs Packmol. Install it with: pixi install -e science"`

Do **not** add a new Python dep for Packmol; shell out to the binary. Keep the
temp-deck writer in a private `_run_packmol(template_atoms, n, box, tol, seed) -> Atoms`.

**Acceptance:**
- `surface_adsorbate` with `element="Cu", facet="fcc111", size=(2,2,3),
  molecule_name="CO", site="ontop"` returns a `Structure` where
  `fragments` has exactly `len-2` entries `== -1` and `2` entries `== 0`,
  `pbc == [True,True,False]`, and `n_fragments == 1`.
- `surface_packing` test **skips if `packmol` binary absent**; when present,
  `n_fragments == cfg.n_molecules` and substrate atoms are all `-1`.

---

## 5. Monte Carlo sampler — `src/traincraft/sampling/monte_carlo.py` (rewrite)

Replace the `NotImplementedError` stub with the full Metropolis sampler. Signature
unchanged: `sample_monte_carlo(structure, calc, job, cfg) -> list[Structure]`.

### 5.1 Setup
- `atoms = structure.atoms.copy(); atoms.calc = calc`
- `rng = np.random.default_rng(cfg.seed)` (fall back to a module default if None;
  do **not** touch the global `np.random`).
- Read fragments via `get_fragments(atoms)`. If `None`: treat all atoms whose
  position is not pinned as a single fragment id `0`, set that array, and
  `logger.warning("monte_carlo: no tc_fragment array; treating system as one rigid fragment")`.
- `mobile = fragment_ids(atoms)`. If empty → raise `ValueError`
  (`"monte_carlo: no mobile fragments to move"`).
- Normalize move probabilities to sum 1.

### 5.2 Energy
- `e_current = atoms.get_potential_energy()` once at start.
- Metropolis: accept Δ if `Δ <= 0` or `rng.random() < exp(-Δ / (kB*T))`,
  `kB = ase.units.kB`, `T = cfg.temperature`.

### 5.3 Moves (operate on a trial copy; only commit on accept)
- **translate**: pick a random mobile fragment, add a uniform vector in
  `[-max_translate, max_translate]^3` to its atoms' positions.
- **rotate**: pick a random mobile fragment, rotate its atoms about their centroid
  by a random axis and an angle uniform in `[-max_rotate, max_rotate]` deg
  (use `scipy.spatial.transform.Rotation` or an explicit Rodrigues rotation).
- **conformer**: only valid if `provenance.extra` (or builder) recorded a SMILES for
  that fragment. Rebuild the fragment with a fresh ETKDG conformer, align its
  centroid + principal axes to the current fragment placement, and substitute the
  coordinates. If RDKit missing or no SMILES recorded → this move is skipped and its
  probability redistributed to translate/rotate (log once). See §5.4.

### 5.4 Conformer-move data dependency
The conformer move needs to know *which SMILES produced each fragment*. Record this
at build time: builders/sources that use a SMILES write
`provenance.extra["fragment_smiles"] = {fragment_id: smiles}` (string keys after
JSON). The MC sampler reads it; absent → conformer move disabled with a warning.
For multi-fragment packings, every fragment shares the same SMILES → map each id to it.

### 5.5 Recording frames
- Every `cfg.interval` steps, snapshot the current accepted `atoms` into a
  `Structure` with `origin="ml_sampled"`, `source="sampler:monte_carlo"`,
  `calculator=<name>`, `parents=[structure.hash]`, and **carry `tc_fragment`
  forward** (it rides along in `atoms.arrays` automatically via `from_ase`).
- If `cfg.refresh_fragments`: after each accepted step, re-run
  `infer_fragments(atoms, scale=cfg.refresh_scale)` and `set_fragments`, so a bond
  forming/breaking re-partitions the system. Log when the fragment count changes.
- Write a trajectory to `job.path("mc.traj")` mirroring `md.py`.
- Return the list of recorded `Structure`s.

**Acceptance:**
- MC on the `surface_adsorbate` Cu+CO system with `calc=EMT`, `steps=50`,
  `interval=10` returns ≥1 frame; **substrate atoms (`-1`) never move** (assert
  their positions identical to the input); each frame carries `tc_fragment`.
- With `refresh_fragments=False` the fragment count is constant across frames.
- Conformer move with no recorded SMILES does not crash (move skipped).
- No `os.chdir`; `os.getcwd()` unchanged after the run (mirror the workspace test).

---

## 6. Wiring & deps

- `geometry/sources/__init__.py`: import `smiles` (guarded — must not hard-fail if
  rdkit missing at *import* time; the registry entry is fine, the ImportError fires
  only when the plugin runs). Follow the rattle pattern: the module imports cleanly,
  the heavy import is inside the function.
- `geometry/builders/__init__.py`: import `surface`.
- `pyproject.toml` + pixi: ensure `science` feature provides `rdkit` (PyPI) and
  `packmol` (conda-forge). Check they are already there; add if not. Do **not** add
  to the `default` env.
- `DESIGN.md` / `ROADMAP.md`: move the implemented items out of "Phase 2/Phase 1
  pending" into done; note fragments as a core concept (§5.1 reference).

---

## 7. Tests (all under `tests/`, must pass in the `science` env; EMT-only ones in `dev`)

| File | Test | Skips if |
|------|------|----------|
| `test_fragments.py` | set/get round-trip; wrong-length rejected; `infer_fragments` on 2 far molecules → 2 ids; framework_mask forces `-1` | — |
| `test_fragments.py` | extxyz round-trip preserves `tc_fragment` (write→read) | — |
| `test_smiles_source.py` | water from `"O"` → 3 atoms, fragments all 0, smiles in provenance.extra | rdkit absent |
| `test_surface_builder.py` | `surface_adsorbate` Cu+CO: fragment counts, pbc, n_fragments==1 | — (EMT/ASE only) |
| `test_surface_builder.py` | `surface_packing` n_fragments==n_molecules; substrate all -1 | packmol binary absent |
| `test_mc_sampler.py` | EMT MC on Cu+CO: ≥1 frame, substrate frozen, fragments carried, cwd unchanged | — |
| `test_config.py` | new models validate; bad adsorbate-count rejected; `builder.file` path resolved relative to config | — |

CI (`dev` env) runs everything that doesn't need rdkit/packmol; the `science` env
job runs the full set. Use `pytest.importorskip("rdkit")` and a
`shutil.which("packmol")` skip guard.

---

## 8. Example configs (under `examples/`)

- `07_co_on_cu_mc.toml` — `surface_adsorbate` (Cu fcc111 + CO) → EMT → `monte_carlo`
  (translate+rotate) → funnel → dataset. **Zero heavy deps** (EMT + ASE), so it runs
  in the default/dev env and doubles as the MC smoke test.
- `08_smiles_molecule.toml` — `smiles = "CCO"` (ethanol) source → EMT MD. Needs
  `science` (rdkit). Header comment says so, mirrors example 02.
- `09_packing_on_surface.toml` — `surface_packing` (Cu + 4×CO2) → EMT → MC with
  `p_conformer > 0`. Needs `science` (packmol + rdkit).

Add `pixi` tasks `example-07/08/09` mirroring the existing ones.

---

## 9. Order of implementation (do in this sequence; commit per numbered group)

1. `core/fragments.py` + `Structure` helpers + `core/__init__` exports + io
   round-trip test. **(foundation — nothing else compiles meaningfully without it)**
2. Config models (§2) + `loader.py` `builder.file` resolution + config tests.
3. `sources/smiles.py` + test.
4. `builders/surface.py` (both builders) + tests.
5. `sampling/monte_carlo.py` rewrite + test.
6. Examples + pixi tasks + DESIGN/ROADMAP doc updates.

Each group is independently testable and reviewable. Keep commit messages to one
sentence per group (per project convention).

---

## 10. Out of scope (explicitly NOT this chunk)

Crystal/defect/slab-only builders (vacancy/substitution/interstitial), layered/2D
stacking, intercalation, polymers (PySoftK), liquids/confined bulk Packmol (only the
surface-coverage case is here), URL/Materials-Project/OPTIMADE/PubChem sources,
strain/rotate/pbc/constraint transforms, SOAP/dscribe diversity. These are later
Phase-1 chunks and must not be started here.
```
