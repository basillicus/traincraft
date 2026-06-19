# Tutorials

The tutorials are step-by-step guides that assume you are new to TrainCraft.
Each one introduces a specific system type or capability and builds on the
previous ones.

You don't have to read them in order — jump to whichever system type is most
relevant to your work.

---

## Learning path

```mermaid
graph LR
    T1["1 · First Dataset<br/>(CNT + EMT + MD)"]
    T2["2 · Molecules<br/>(SMILES + RDKit)"]
    T3["3 · Surfaces<br/>(adsorbate + MC)"]
    T4["4 · Crystals<br/>(defects)"]
    T5["5 · Slabs<br/>(strain + transforms)"]
    T6["6 · 2D Materials<br/>(graphene, MoS₂)"]
    T7["7 · MACE-MP0<br/>(foundation model)"]
    T8["8 · Selection Funnel<br/>(diversity, budget)"]
    T9["9 · Interop<br/>(pymatgen, RDKit)"]

    T1 --> T2
    T1 --> T4
    T2 --> T3
    T4 --> T5
    T4 --> T6
    T1 --> T7
    T1 --> T8
    T3 --> T9
```

---

## Quick reference

| Tutorial | System type | Extra deps | Pixi env |
|---|---|---|---|
| [1 · First Dataset](01-first-dataset.md) | Carbon nanotube | None | `default` |
| [2 · Molecules](02-molecules.md) | Small molecules, SMILES | RDKit | `science` |
| [3 · Surfaces](03-surfaces.md) | Adsorbate + MC sampler | (RDKit optional) | `default` / `science` |
| [4 · Crystals](04-crystals-defects.md) | Bulk + defects | None | `default` |
| [5 · Slabs & Strain](05-slabs-strain.md) | Slab + transforms | None | `default` |
| [6 · 2D Materials](06-2d-materials.md) | Graphene, hBN, MoS₂ | None | `default` |
| [7 · MACE-MP0](07-mace-mp0.md) | Any periodic system | torch, mace-torch | `mace` |
| [8 · Selection Funnel](08-selection-funnel.md) | Any | None | `default` |
| [9 · Interop](09-interop.md) | Any | pymatgen, RDKit | `science` |
