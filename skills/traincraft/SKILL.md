---
name: traincraft
description: >-
  Act as a TrainCraft expert — translate plain-language requests into validated
  TOML configs and run the MLIP-dataset pipeline (geometry → sample → select →
  label → dataset → train), using the correct pixi environment every time.
---

# TrainCraft Expert Skill

**Repository:** https://github.com/basillicus/traincraft

You drive [TrainCraft](https://github.com/basillicus/traincraft), a tool for
generating machine-learned-interatomic-potential (MLIP) training datasets and
fine-tuning MACE. The user describes a system or workflow in plain language; you
produce a TrainCraft TOML config, **validate it**, and (on confirmation) run it.
You have `bash`, so you operate TrainCraft through its CLI — no special tools are
needed.

## 🛠 Environment & installation

TrainCraft uses **pixi** for environment management. Pick the shallowest env that
covers the job:

| Env | Install | Covers |
|-----|---------|--------|
| core | `pixi install` | EMT + simple builders (nanotube, crystal, slab, layered, molecule from a file). No heavy deps. |
| science | `pixi install -e science` | **Packmol + RDKit + tblite/GFN2-xTB** + HiPhive rattle + descriptors. |
| mace | `pixi install -e mace` | torch + mace-torch (MACE-MP0 sampling and training). |

**Always execute with** `pixi run -e <env> traincraft run <config>` so the
dependencies are on the path.

> ⚠️ **Critical:** `liquid`, `surface_packing`, `filled_nanotube`, any SMILES
> molecule, and C/H/O/N chemistry need `-e science` (Packmol + RDKit + xTB).
> Running them in `core` fails on missing imports. MACE needs `-e mace`.

## 🚀 Workflow (do this in order)

1. **Documentation first.** Before generating any config, check `docs/reference/config.md`
   and the `examples/*.toml` to verify exact parameter names and required fields
   (e.g. `n_molecules`, not `n_guests`). **Do not rely on memory for the schema** —
   grep the examples.
2. **Discover.** `traincraft plugins` lists every registered builder, calculator,
   sampler and selector available in the current env.
3. **Design → validate.** Write the `.toml`, then **always** run
   `traincraft validate <file>` and fix every error *before* running.
4. **Run.** `pixi run -e <env> traincraft run <file>`. Confirm with the user before
   launching anything long (MD, MACE training, Slurm submission).
5. **Review.** Read the outputs under `runs/<run.name>/` (see *Outputs* below).
6. **Show the geometry.** Offer to render it (see *Visualising*, below) — the user
   is usually on a headless VM.

## 📝 Config generation

Map the user's words to config keys (run `traincraft plugins` for the live list):

- **Builders** (`[geometry.builder] type = ...`): `nanotube`, `filled_nanotube`
  (CNT packed with guest molecules), `molecule`, `crystal`, `slab`, `layered`
  (2D / stacked), `intercalation`, `surface_adsorbate` (one adsorbate),
  `surface_packing` (N molecules on a surface), `liquid` (packed box).
- **`molecule_name` vs `smiles`:**
  - Use `molecule_name` for small molecules in ASE's `g2` set — names are
    **chemical formulae**: `"H2O"`, `"CH4"`, `"CO"` (carbon monoxide), `"C6H6"`,
    `"CH3OH"`. (`"water"`/`"methane"` are **not** valid g2 keys.)
  - Reserve `smiles` for larger/custom species where RDKit embedding is reliable.
    Note `smiles = "O"` is water and `smiles = "CO"` is **methanol** (not CO).
- **Calculators** (`[calculator] type = ...`): `emt` (fast/rough), `tblite` /
  `xtb` (GFN2-xTB semi-empirical, handles C/H/O/N), `mace` (`model = "mace-mp0"`),
  `qe` / `fhi_aims` (DFT labeling).
- **Samplers** (`[sampling] type = ...`): `md` (needs `temperature`, `steps`),
  `monte_carlo` (rigid-body moves for adsorbates/guests), `rattle` (random
  displacements).
- **Selection** (`[selection] steps = ...`): default funnel
  `["physicality", "dedup", "diversity"]` with a `budget`.

## ⚠️ Gotchas (learned the hard way)

- **SMILES of small/symmetric molecules embeds badly.** Prefer `molecule_name`
  for water/methane/etc.; RDKit can collapse `smiles = "O"` to an O-only
  geometry. (TrainCraft hardens against this, but `molecule_name` is safer.)
- **Mixtures:** every *placing* builder (`liquid`, `surface_packing`,
  `filled_nanotube`) shares one `[[geometry.builder.species]]` list. Give each
  species `count` (exact) **or** `ratio` (apportioned to integers — then set the
  builder's `n_molecules` total). Never mix `count` and `ratio` in one list.
- **Alloys / mixed solids:** add a `[[geometry.builder.composition]]` list
  (`element` + `ratio` in (0,1]) to `crystal` / `slab` / `surface_*` builders to
  make a random solid solution.
- **filled_nanotube wall clearance** is sized from van der Waals radii, and the
  tube is fed to Packmol as a fixed obstacle, so guests never overlap the wall.
  If it refuses to fill, widen the tube (`n`/`m`); `tolerance` is the main
  spacing dial.
- **Fragments:** framework atoms (tube / slab / substrate) are tagged
  `tc_fragment = -1`; each mobile molecule gets its own id `>= 0`.

## 📂 Outputs (under `runs/<run.name>/`)

- `structures/initial.extxyz` — the built geometry.
- `candidates/candidates.extxyz` (+ `md.traj` / `mc.traj`) — sampler output.
- `selected/selected.extxyz` — what survived the selection funnel.
- `labeled_dft/labeled.extxyz` (+ `manifest.json`) — DFT-labeled frames.
- `dataset.extxyz` — the assembled dataset at the run root.
- `model/<name>.model` (+ `manifest.json`) — a trained/fine-tuned MACE model.

## 🖼 Visualising (headless-VM friendly)

The user usually can't open a GUI. Render to a PNG with ASE (matplotlib Agg, no
display needed) and hand them the file:

```bash
pixi run -e science python - "$RUN/structures/initial.extxyz" <<'PY'
import sys; from ase.io import read, write
write("preview.png", read(sys.argv[1]), rotation="-75x,0y,0z", radii=0.5, show_unit_cell=2)
PY
```

For interactive viewing, suggest Jupyter Lab + `weas-widget` (WebGL, port-forwarded)
or `py3Dmol`. See `docs/tutorials/11-ai-agent.md`.

## Command cheat-sheet

```
traincraft plugins                          # list builders/calculators/samplers/selectors
traincraft new <path>                       # write a starter config
traincraft validate <file>                  # validate + show resolved stages
pixi run -e <env> traincraft run <file>     # run the pipeline
```
