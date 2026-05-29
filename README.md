# TrainCraft

Modular generation of training datasets for — and training of — machine-learned
interatomic potentials (MLIPs), with a MACE-centric **active-learning** loop that
only spends expensive DFT on structures that actually improve the model.

> **Status:** greenfield rewrite in progress. The architecture is in
> [`DESIGN.md`](DESIGN.md) and the phased plan in [`ROADMAP.md`](ROADMAP.md).
> The legacy code (used for the Raman work) is preserved at git tag/commit
> `e39647c`. This rewrite has **no backward compatibility** with it.

## What works today (Phase 0 + vertical slice)

A complete, if thin, pipeline you can run with **zero heavy dependencies**:

```
geometry builder → cheap calculator → sampler → selection funnel → dataset
```

```bash
pip install -e ".[dev]"
traincraft run examples/walking_skeleton.toml
```

This builds a small carbon nanotube, runs a short MD with ASE's EMT calculator,
filters the frames through the selection funnel (physicality → dedup →
diversity), and writes a provenance-tagged `dataset.extxyz` under
`runs/skeleton_cnt/`.

### CLI

```bash
traincraft run <config.toml>       # run the whole workflow
traincraft validate <config.toml>  # validate and show resolved stages
traincraft new <path.toml>         # scaffold a starter config
traincraft plugins                 # list registered plugins
```

## One TOML drives everything

A single validated TOML declares the entire workflow; `traincraft run` executes
it end to end. Adding a new geometry builder, calculator, sampler, or selector is
one self-registering file — see [`DESIGN.md`](DESIGN.md).

## Use it as a library

Every stage is a pure function, usable on its own:

```python
import traincraft as tc

cfg = tc.load_config("examples/walking_skeleton.toml")
summary = tc.run_pipeline(cfg)

# or compose pieces directly
structure = tc.build_geometry(cfg.geometry)
calc = tc.make_calculator(cfg.calculator)
```

## Environment (pixi)

Pixi mixes conda-forge (packmol, hiphive, DFT tooling) and PyPI (mace-torch,
pysoftk, …) in one lockfile:

```bash
pixi install            # default env
pixi run test           # pytest
pixi run slice          # run the walking-skeleton example
```

Optional capability groups: `mace`, `geometry`, `semiempirical`, `descriptors`,
`sampling` (HiPhive), `dev`.

## License

See [`LICENSE`](LICENSE). If you use the legacy TrainCraft in research, cite the
Zenodo DOI `10.5281/zenodo.8174842` and the packages it builds on.
