# Write a Custom Builder

Adding a new geometry builder to TrainCraft takes three steps: write the
function, add a config model, and register both. No existing files need to
be modified beyond the config union and the `__init__` import.

---

## Step 1 — Write the builder function

Create `src/traincraft/geometry/builders/my_builder.py`:

```python
from __future__ import annotations

from ...core import Provenance, Structure, register


@register("builder", "my_builder")   # (1)
def build_my_builder(cfg) -> Structure:   # (2)
    # Build an ase.Atoms object using cfg fields
    from ase import Atoms
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, cfg.bond_length]])
    atoms.center(vacuum=cfg.vacuum / 2)

    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated",
            source=f"builder:my_builder:H2-{cfg.bond_length}",
        ),
    )
```

1. **`@register("builder", "my_builder")`** — registers the function under the
   key `"my_builder"`. This is the string you use in `type = "my_builder"` in TOML.

2. **`cfg`** — the config model instance. It has whatever fields you declare in
   the pydantic model (Step 2).

---

## Step 2 — Add a config model

Add to `src/traincraft/config/models.py`:

```python
class MyBuilder(TCModel):
    type: Literal["my_builder"] = "my_builder"
    bond_length: float = 0.74      # Å
    vacuum: float = 6.0

    @model_validator(mode="after")
    def _check(self) -> MyBuilder:
        if self.bond_length <= 0:
            raise ValueError("bond_length must be positive")
        return self
```

Then add `MyBuilder` to the `BuilderConfig` union:

```python
BuilderConfig = Annotated[
    NanotubeBuilder | MoleculeBuilder | ... | MyBuilder,   # ← add here
    Field(discriminator="type"),
]
```

---

## Step 3 — Register on import

Add to `src/traincraft/geometry/builders/__init__.py`:

```python
from . import (
    crystal, layered, molecule, my_builder, nanotube, slab, surface,  # ← add
)
```

---

## Step 4 — Use it in a TOML

```toml
[geometry.builder]
type        = "my_builder"
bond_length = 0.80
vacuum      = 8.0
```

---

## Step 5 — Write a test

```python
from traincraft.config.models import GeometryConfig
from traincraft.geometry import build_geometry
from traincraft.geometry.builders.my_builder import MyBuilder

def test_my_builder_atom_count():
    s = build_geometry(GeometryConfig(builder=MyBuilder(bond_length=0.74)))
    assert len(s.atoms) == 2
    assert s.provenance.source.startswith("builder:my_builder")
```

---

## Tips

!!! tip "Fragment tagging"
    If your builder creates a system with distinct mobile components (adsorbates,
    molecules, etc.), set the `tc_fragment` array so the MC sampler can use it:

    ```python
    import numpy as np
    from ...core import set_fragments, FRAMEWORK

    frag = np.full(len(atoms), FRAMEWORK, dtype=int)
    frag[n_substrate:] = 0   # mobile fragment 0
    set_fragments(atoms, frag)
    ```

!!! tip "External tool integration"
    If your builder calls an external binary, check that it's on `PATH` and
    raise `ImportError` with an install hint if not — matching the pattern in
    `surface.py` (`shutil.which("packmol")`).

!!! tip "Config validation"
    Use `@model_validator` to catch logical errors early (e.g., mutually exclusive
    options). The `TCModel` base class already sets `extra="forbid"`, so unknown
    fields fail immediately without any extra work.
