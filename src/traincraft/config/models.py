"""Validated, typed configuration models (pydantic v2).

A single TOML file declares the whole workflow. Plugin sections are
discriminated unions on a ``type`` field, so adding a plugin = adding a model
here plus a registry entry. ``extra="forbid"`` makes typos fail fast.
"""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class TCModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


# --------------------------------------------------------------------------- run
class RunConfig(TCModel):
    name: str = "traincraft_run"
    outdir: str = "runs"
    seed: int | None = 42


# ------------------------------------------------------------------------ sources
class FileSource(TCModel):
    type: Literal["file"] = "file"
    path: str  # any ASE-readable format


class ScratchSource(TCModel):
    type: Literal["scratch"] = "scratch"
    molecule: str | None = None  # ase.build.molecule name (e.g. "H2O")
    bulk: str | None = None  # ase.build.bulk symbol (e.g. "Cu")


SourceConfig = Annotated[FileSource | ScratchSource, Field(discriminator="type")]


# ----------------------------------------------------------------------- builders
class NanotubeBuilder(TCModel):
    type: Literal["nanotube"] = "nanotube"
    n: int = 8
    m: int = 0
    length: int = 1
    bond: float = 1.42
    vacuum: float = 6.0


class MoleculeBuilder(TCModel):
    type: Literal["molecule"] = "molecule"
    name: str | None = None  # ASE g2 name (e.g. "CO2")
    smiles: str | None = None  # RDKit path (geometry chunk)
    vacuum: float = 6.0

    @model_validator(mode="after")
    def _need_one(self) -> MoleculeBuilder:
        if not self.name and not self.smiles:
            raise ValueError("molecule builder needs 'name' or 'smiles'")
        return self


BuilderConfig = Annotated[
    NanotubeBuilder | MoleculeBuilder, Field(discriminator="type")
]


# --------------------------------------------------------------------- transforms
class VacuumTransform(TCModel):
    type: Literal["vacuum"] = "vacuum"
    amount: float = 6.0


class SupercellTransform(TCModel):
    type: Literal["supercell"] = "supercell"
    repeat: tuple[int, int, int] = (1, 1, 1)


class PerturbTransform(TCModel):
    type: Literal["perturb"] = "perturb"
    stddev: float = 0.05


TransformConfig = Annotated[
    VacuumTransform | SupercellTransform | PerturbTransform,
    Field(discriminator="type"),
]


# -------------------------------------------------------------------------- geometry
class GeometryConfig(TCModel):
    source: SourceConfig | None = None
    builder: BuilderConfig | None = None
    transforms: list[TransformConfig] = Field(default_factory=list)

    @model_validator(mode="after")
    def _need_one(self) -> GeometryConfig:
        if self.source is None and self.builder is None:
            raise ValueError("geometry needs a 'source' or a 'builder'")
        return self


# ------------------------------------------------------------------------ calculators
class EmtCalc(TCModel):
    type: Literal["emt"] = "emt"


class TbliteCalc(TCModel):
    type: Literal["tblite"] = "tblite"
    method: str = "GFN2-xTB"


class XtbCalc(TCModel):
    type: Literal["xtb"] = "xtb"
    method: str = "GFN2-xTB"


class MaceCalc(TCModel):
    type: Literal["mace"] = "mace"
    model: str = "mace-mp0"  # mace-mp0 | mace-off23
    model_path: str | None = None  # local fine-tuned model
    device: str = "cpu"
    default_dtype: str = "float32"


CalculatorConfig = Annotated[
    EmtCalc | TbliteCalc | XtbCalc | MaceCalc, Field(discriminator="type")
]


# -------------------------------------------------------------------------- sampling
class MdSampling(TCModel):
    type: Literal["md"] = "md"
    temperature: float = 500.0
    steps: int = 1000
    interval: int = 20
    timestep: float = 1.0  # fs
    friction: float = 0.02


class RattleSampling(TCModel):
    type: Literal["rattle"] = "rattle"
    method: Literal["mc", "standard"] = "mc"
    n_structures: int = 10
    std: float = 0.12
    min_distance: float = 1.3


class MonteCarloSampling(TCModel):
    type: Literal["monte_carlo"] = "monte_carlo"
    steps: int = 1000
    temperature: float = 500.0
    conformers: bool = False  # RDKit conformer moves (next chunk)


SamplingConfig = Annotated[
    MdSampling | RattleSampling | MonteCarloSampling, Field(discriminator="type")
]


# ------------------------------------------------------------------------- selection
class SelectionConfig(TCModel):
    steps: list[str] = Field(default_factory=lambda: ["physicality", "dedup", "diversity"])
    budget: int | None = 50
    min_distance: float = 0.7  # physicality: min interatomic distance (Angstrom)


# --------------------------------------------------------------------------- dataset
class DatasetConfig(TCModel):
    path: str = "dataset"
    format: Literal["extxyz"] = "extxyz"


# ------------------------------------------------------------------------------ root
class TrainCraftConfig(TCModel):
    run: RunConfig = Field(default_factory=RunConfig)
    geometry: GeometryConfig | None = None
    calculator: CalculatorConfig | None = None
    sampling: SamplingConfig | None = None
    selection: SelectionConfig | None = None
    dataset: DatasetConfig | None = None
