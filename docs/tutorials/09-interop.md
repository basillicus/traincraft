# Tutorial 9 · ASE · pymatgen · RDKit

**What you'll learn:** how to convert TrainCraft structures to and from pymatgen
and RDKit — enabling you to use materials databases, symmetry analysis,
cheminformatics, and external tools alongside TrainCraft.

**Prerequisites:** Tutorial 1. Requires `pixi install -e science`.

---

## The converter module

`traincraft.core.converter` provides four functions:

```python
from traincraft.core.converter import (
    ase_to_pymatgen,   # ase.Atoms → pymatgen Structure or Molecule
    pymatgen_to_ase,   # pymatgen Structure/Molecule → ase.Atoms
    ase_to_rdkit,      # ase.Atoms → rdkit.Chem.Mol (non-periodic only)
    rdkit_to_ase,      # rdkit.Chem.Mol → ase.Atoms
)
```

These are also available as methods on the `Structure` class:

```python
from traincraft import Structure

s = Structure.from_ase(atoms)
pmg  = s.to_pymatgen()    # → pymatgen Structure or Molecule
rdmol = s.to_rdkit()      # → rdkit.Chem.Mol (non-periodic only)

s2 = Structure.from_pymatgen(pmg)
s3 = Structure.from_rdkit(rdmol)
```

---

## ASE ↔ pymatgen

### Periodic structure (crystal/slab)

```python
from ase.build import bulk
from traincraft.core.converter import ase_to_pymatgen, pymatgen_to_ase

cu = bulk("Cu", "fcc", a=3.61, cubic=True)
pmg = ase_to_pymatgen(cu)   # → pymatgen.core.Structure

print(type(pmg).__name__)    # "Structure"
print(pmg.composition)       # Comp: Cu4
print(pmg.lattice)           # Lattice(a=3.61, b=3.61, c=3.61, ...)

# Round-trip
cu_back = pymatgen_to_ase(pmg)
```

### Non-periodic molecule

```python
from ase.build import molecule
from traincraft.core.converter import ase_to_pymatgen

h2o = molecule("H2O")
mol = ase_to_pymatgen(h2o)   # → pymatgen.core.Molecule (not Structure)
print(type(mol).__name__)     # "Molecule"
print(mol.formula)            # "H2 O1"
```

!!! note "Partially periodic (slab)"
    A slab periodic in xy but not z is treated as **non-periodic** by pymatgen
    (it becomes a `Molecule`). pymatgen has no direct analogue for partial
    periodicity. If you need a pymatgen `Structure` for a slab, set `pbc=True`
    first (using `set_pbc` transform), perform your analysis, then reset.

---

## Use case: querying the Materials Project

pymatgen provides a client for the Materials Project API. Once you fetch a
structure from MP, bring it into TrainCraft with `Structure.from_pymatgen`:

```python
from mp_api.client import MPRester
from traincraft import Structure
from traincraft.core.provenance import Provenance

with MPRester("YOUR_API_KEY") as mpr:
    docs = mpr.materials.summary.search(
        material_ids=["mp-30"],   # Cu (mp-30)
        fields=["structure"],
    )
    pmg_struct = docs[0].structure

tc_struct = Structure.from_pymatgen(
    pmg_struct,
    provenance=Provenance(
        origin="generated",
        source="mp:mp-30",
        extra={"mp_id": "mp-30"},
    ),
)
print(tc_struct.atoms.get_chemical_symbols()[:4])
```

---

## ASE ↔ RDKit

### Convert to RDKit Mol

RDKit conversion is only supported for **non-periodic** structures. The
converter writes the atoms to an in-memory XYZ block, parses it with
`Chem.MolFromXYZBlock`, and then uses `DetermineBonds` (the xyz2mol algorithm)
to perceive bonding.

```python
from ase.build import molecule
from traincraft.core.converter import ase_to_rdkit

h2o = molecule("H2O")
mol = ase_to_rdkit(h2o)     # charge=0 by default

print(mol.GetNumAtoms())    # 3
print(mol.GetNumBonds())    # 2

# RDKit SMILES
from rdkit.Chem import MolToSmiles
print(MolToSmiles(mol))     # "O"
```

If your molecule has a non-zero charge:

```python
mol = ase_to_rdkit(cation_atoms, charge=+1)
```

### Convert from RDKit Mol

If you have an RDKit `Mol` with embedded 3D coordinates, convert it back:

```python
from rdkit.Chem import MolFromSmiles, AddHs
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule, ETKDGv3
from traincraft.core.converter import rdkit_to_ase

mol = AddHs(MolFromSmiles("CCO"))     # ethanol
params = ETKDGv3()
params.randomSeed = 42
EmbedMolecule(mol, params)
MMFFOptimizeMolecule(mol)

atoms = rdkit_to_ase(mol)             # first conformer by default
# atoms = rdkit_to_ase(mol, conf_id=1)  # specific conformer
```

---

## Use case: SMILES → TrainCraft → MACE

A complete pipeline: start from a SMILES string, generate conformers with RDKit,
build a TrainCraft structure, and run MACE-OFF23:

```python
from rdkit import Chem
from rdkit.Chem import AllChem, AddHs
from traincraft import Structure
from traincraft.core.provenance import Provenance
from traincraft.core.converter import rdkit_to_ase
from traincraft.config.models import MaceCalc

# Generate 5 conformers
smiles = "c1ccccc1"   # benzene
mol = AddHs(Chem.MolFromSmiles(smiles))
params = AllChem.ETKDGv3()
AllChem.EmbedMultipleConfs(mol, numConfs=5, params=params)
AllChem.MMFFOptimizeMoleculeConfs(mol)

structures = [
    Structure.from_rdkit(mol, conf_id=i,
                         provenance=Provenance(origin="generated",
                                               source=f"smiles:{smiles}:conf{i}"))
    for i in range(mol.GetNumConformers())
]
print(f"Generated {len(structures)} conformers")
```

---

## Use case: symmetry analysis with pymatgen

```python
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from traincraft.geometry import build_geometry
from traincraft.config.models import CrystalBuilder, GeometryConfig

s = build_geometry(GeometryConfig(builder=CrystalBuilder(
    name="Si", crystalstructure="diamond", a=5.43, cubic=True
)))

pmg = s.to_pymatgen()
sga = SpacegroupAnalyzer(pmg)
print(sga.get_space_group_symbol())       # "Fd-3m"
print(sga.get_space_group_number())       # 227
primitive = sga.get_primitive_standard_structure()
print(len(primitive))                     # 2 atoms
```

---

## Summary

| Conversion | Function | Notes |
|---|---|---|
| ASE → pymatgen | `ase_to_pymatgen` | Periodic → `Structure`; molecular → `Molecule` |
| pymatgen → ASE | `pymatgen_to_ase` | Works for both `Structure` and `Molecule` |
| ASE → RDKit | `ase_to_rdkit` | Non-periodic only; bonds perceived by xyz2mol |
| RDKit → ASE | `rdkit_to_ase` | Uses `conf_id` to select conformer |

**Next:** [Tutorial 10](10-training.md) — train a MACE model on the dataset you've
built. Or check out the [Concepts](../concepts/architecture.md) section for a
deeper understanding of how TrainCraft works under the hood, or the
[Config Schema](../reference/config.md) for the complete field reference.
