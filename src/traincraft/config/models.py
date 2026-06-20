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


class SmilesSource(TCModel):
    type: Literal["smiles"] = "smiles"
    smiles: str
    n_conformers: int = 1
    optimize: bool = True
    seed: int | None = None
    vacuum: float = 6.0


class UrlSource(TCModel):
    type: Literal["url"] = "url"
    url: str
    format: str | None = None  # ASE format override; inferred from suffix if None
    timeout: float = 30.0


class MaterialsProjectSource(TCModel):
    type: Literal["materials_project"] = "materials_project"
    material_id: str  # e.g. "mp-149"
    api_key: str | None = None  # else taken from $MP_API_KEY
    conventional: bool = True  # conventional (vs primitive) unit cell


class OptimadeSource(TCModel):
    type: Literal["optimade"] = "optimade"
    base_url: str  # an OPTIMADE provider base URL
    filter: str | None = None  # OPTIMADE filter; first hit is returned
    timeout: float = 30.0


class PubchemSource(TCModel):
    type: Literal["pubchem"] = "pubchem"
    # exactly one identifier; resolved to a 3D conformer via PubChem PUG REST
    cid: int | None = None
    name: str | None = None
    smiles: str | None = None
    timeout: float = 30.0

    @model_validator(mode="after")
    def _one_id(self) -> PubchemSource:
        if sum(x is not None for x in (self.cid, self.name, self.smiles)) != 1:
            raise ValueError("pubchem source needs exactly one of 'cid', 'name', or 'smiles'")
        return self


SourceConfig = Annotated[
    FileSource
    | ScratchSource
    | SmilesSource
    | UrlSource
    | MaterialsProjectSource
    | OptimadeSource
    | PubchemSource,
    Field(discriminator="type"),
]


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


_FACETS = Literal[
    "fcc111", "fcc100", "fcc110", "bcc110", "bcc100", "bcc111", "hcp0001"
]


class DefectSpec(TCModel):
    kind: Literal["vacancy", "substitution", "interstitial"]
    index: int | None = None  # target atom (vacancy/substitution); default 0
    element: str | None = None  # new element (substitution/interstitial)
    position: tuple[float, float, float] | None = None  # interstitial site
    cartesian: bool = False  # interstitial position is Cartesian (else fractional)

    @model_validator(mode="after")
    def _check(self) -> DefectSpec:
        if self.kind == "substitution" and self.element is None:
            raise ValueError("substitution defect needs 'element'")
        if self.kind == "interstitial" and (
            self.element is None or self.position is None
        ):
            raise ValueError("interstitial defect needs 'element' and 'position'")
        return self


class CrystalBuilder(TCModel):
    type: Literal["crystal"] = "crystal"
    name: str  # ase.build.bulk name, e.g. "Cu", "NaCl", "Si"
    crystalstructure: str | None = None  # fcc/bcc/hcp/diamond/rocksalt/...
    a: float | None = None
    b: float | None = None
    c: float | None = None
    cubic: bool = False
    orthorhombic: bool = False
    supercell: tuple[int, int, int] = (1, 1, 1)
    defects: list[DefectSpec] = Field(default_factory=list)


class SlabBuilder(TCModel):
    type: Literal["slab"] = "slab"
    element: str
    # Mode A: named facet (element-only). Mode B: miller indices from a bulk.
    facet: _FACETS | None = "fcc111"
    miller: tuple[int, int, int] | None = None
    crystalstructure: str | None = None  # for miller mode (fcc/bcc/hcp/...)
    a: float | None = None
    c: float | None = None
    cubic: bool = False
    layers: int = 4
    size: tuple[int, int, int] = (3, 3, 4)  # (nx, ny, layers); layers used by facet mode
    vacuum: float = 12.0
    orthogonal: bool = False  # facet mode
    periodic: bool = False  # miller mode: keep periodic along the surface normal

    @model_validator(mode="after")
    def _one_mode(self) -> SlabBuilder:
        if self.miller is not None:
            self.facet = None  # miller mode takes precedence; clear the default facet
        elif self.facet is None:
            raise ValueError("slab builder needs either 'facet' or 'miller'")
        return self


class SurfaceAdsorbateBuilder(TCModel):
    type: Literal["surface_adsorbate"] = "surface_adsorbate"
    # substrate
    element: str
    facet: _FACETS = "fcc111"
    size: tuple[int, int, int] = (3, 3, 4)
    vacuum: float = 12.0
    # adsorbate: exactly one of molecule_name | smiles | file
    molecule_name: str | None = None  # ase.build.molecule g2 name
    smiles: str | None = None
    file: str | None = None  # path resolved relative to config file
    # placement
    site: Literal["ontop", "bridge", "hollow", "fcc", "hcp"] = "ontop"
    height: float = 2.0
    offset: tuple[float, float] | None = None

    @model_validator(mode="after")
    def _one_adsorbate(self) -> SurfaceAdsorbateBuilder:
        n = sum(x is not None for x in (self.molecule_name, self.smiles, self.file))
        if n != 1:
            raise ValueError(
                "surface_adsorbate needs exactly one of 'molecule_name', 'smiles', or 'file'"
            )
        return self


class SurfacePackingBuilder(TCModel):
    type: Literal["surface_packing"] = "surface_packing"
    element: str
    facet: _FACETS = "fcc111"
    size: tuple[int, int, int] = (4, 4, 4)
    vacuum: float = 20.0
    # adsorbate molecules to pack above the slab: exactly one
    molecule_name: str | None = None
    smiles: str | None = None
    file: str | None = None
    # packing parameters
    n_molecules: int = 4
    tolerance: float = 2.0
    region_height: float = 8.0
    gap: float = 2.0
    seed: int | None = None

    @model_validator(mode="after")
    def _one_adsorbate(self) -> SurfacePackingBuilder:
        n = sum(x is not None for x in (self.molecule_name, self.smiles, self.file))
        if n != 1:
            raise ValueError(
                "surface_packing needs exactly one of 'molecule_name', 'smiles', or 'file'"
            )
        return self


class LayeredBuilder(TCModel):
    type: Literal["layered"] = "layered"
    material: Literal["graphene", "hbn", "mx2"] = "graphene"
    formula: str | None = None  # mx2 only, e.g. "MoS2", "WSe2"
    a: float | None = None  # in-plane lattice constant override
    size: tuple[int, int] = (1, 1)  # in-plane repeats of the monolayer
    n_layers: int = 2
    interlayer_spacing: float = 3.35  # Å (graphite default)
    stacking: Literal["AA", "AB"] = "AB"
    twist: float = 0.0  # degrees; nonzero -> non-periodic moiré flake
    vacuum: float = 15.0


class LiquidSpecies(TCModel):
    # exactly one of molecule_name | smiles | file identifies the species
    molecule_name: str | None = None
    smiles: str | None = None
    file: str | None = None
    count: int = 1

    @model_validator(mode="after")
    def _one_source(self) -> LiquidSpecies:
        if sum(x is not None for x in (self.molecule_name, self.smiles, self.file)) != 1:
            raise ValueError(
                "liquid species needs exactly one of 'molecule_name', 'smiles', 'file'"
            )
        return self


class LiquidBuilder(TCModel):
    type: Literal["liquid"] = "liquid"
    species: list[LiquidSpecies]
    box: tuple[float, float, float] | None = None  # explicit cell edges (Å)
    density: float | None = None  # g/cm³; cubic box derived when 'box' is unset
    tolerance: float = 2.0  # Packmol min separation / wall inset (Å)
    pbc: bool = True
    seed: int | None = None

    @model_validator(mode="after")
    def _box_or_density(self) -> LiquidBuilder:
        if not self.species:
            raise ValueError("liquid builder needs at least one species")
        if (self.box is None) == (self.density is None):
            raise ValueError("liquid builder needs exactly one of 'box' or 'density'")
        return self


class IntercalationBuilder(TCModel):
    type: Literal["intercalation"] = "intercalation"
    host: LayeredBuilder  # planar layered host (graphene/hbn); mx2 rejected
    guest: str = "Li"  # intercalant element symbol
    n_per_gallery: int = 1  # guest atoms per filled gallery (in-plane grid)
    stage: int = 1  # fill galleries whose index % stage == 0
    gallery_expansion: float = 0.0  # Å added to c per filled gallery


BuilderConfig = Annotated[
    NanotubeBuilder
    | MoleculeBuilder
    | SurfaceAdsorbateBuilder
    | SurfacePackingBuilder
    | CrystalBuilder
    | SlabBuilder
    | LayeredBuilder
    | LiquidBuilder
    | IntercalationBuilder,
    Field(discriminator="type"),
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


class StrainTransform(TCModel):
    type: Literal["strain"] = "strain"
    # exactly one of: hydrostatic OR voigt (e_xx, e_yy, e_zz, e_yz, e_xz, e_xy)
    hydrostatic: float | None = None
    voigt: tuple[float, float, float, float, float, float] | None = None

    @model_validator(mode="after")
    def _one_form(self) -> StrainTransform:
        n = (self.hydrostatic is not None) + (self.voigt is not None)
        if n != 1:
            raise ValueError("strain needs exactly one of 'hydrostatic' or 'voigt'")
        return self


class RotateTransform(TCModel):
    type: Literal["rotate"] = "rotate"
    angle: float  # degrees
    axis: Literal["x", "y", "z"] | tuple[float, float, float] = "z"
    rotate_cell: bool = False


class SetPbcTransform(TCModel):
    type: Literal["set_pbc"] = "set_pbc"
    pbc: bool | tuple[bool, bool, bool] = True


class ConstraintsTransform(TCModel):
    type: Literal["constraints"] = "constraints"
    # selectors are OR-ed; at least one must select an atom
    indices: list[int] | None = None
    elements: list[str] | None = None  # fix all atoms of these symbols
    fragments: list[int] | None = None  # fix atoms with these tc_fragment ids
    below_z: float | None = None  # fix atoms with Cartesian z below this (Å)


TransformConfig = Annotated[
    VacuumTransform
    | SupercellTransform
    | PerturbTransform
    | StrainTransform
    | RotateTransform
    | SetPbcTransform
    | ConstraintsTransform,
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


# Properties a DFT labeler may be asked to produce beyond E/F/stress (always on).
_DFT_PROP = Literal["dipole", "polarizability"]


class FhiAimsCalc(TCModel):
    """FHI-aims DFT labeler (ASE ``FileIOCalculator``).

    The run command is *not* configured here: it is injected from the
    environment (``TRAINCRAFT_AIMS_COMMAND``, default ``aims.x``) so the plugin
    stays container-agnostic (see DESIGN §20.3). Likewise the species directory
    falls back to ``TRAINCRAFT_AIMS_SPECIES_DIR`` / ``AIMS_SPECIES_DIR`` when
    ``species_dir`` is left unset.
    """

    type: Literal["fhi_aims"] = "fhi_aims"
    # extra properties to label (energy/forces/stress are always requested)
    properties: list[_DFT_PROP] = Field(default_factory=list)
    # core DFT settings
    xc: str = "pbe"
    species_defaults: str = "light"  # basis level: light | intermediate | tight | ...
    species_dir: str | None = None  # explicit path; else taken from env
    kpts: tuple[int, int, int] | None = None  # Monkhorst-Pack grid (periodic only)
    relativistic: str = "atomic_zora scalar"
    spin: Literal["none", "collinear"] = "none"
    # periodicity hint for DFPT mode selection; auto-True when kpts is given
    periodic: bool | None = None
    # arbitrary control.in keywords passed through verbatim to ASE-Aims
    extra: dict = Field(default_factory=dict)


class QeCalc(TCModel):
    """Quantum ESPRESSO ``pw.x`` DFT labeler (ASE ``Espresso``).

    Command injected from ``TRAINCRAFT_PW_COMMAND`` (default ``pw.x``);
    pseudopotential directory from ``pseudo_dir`` or
    ``TRAINCRAFT_PW_PSEUDO_DIR`` / ``ESPRESSO_PSEUDO``.
    """

    type: Literal["qe"] = "qe"
    # only "dipole" is supported via SCF; polarizability needs ph.x (see dft.py)
    properties: list[_DFT_PROP] = Field(default_factory=list)
    pseudo_dir: str | None = None  # else taken from env
    pseudopotentials: dict = Field(default_factory=dict)  # {symbol: UPF filename}
    kpts: tuple[int, int, int] | None = None
    kspacing: float | None = None  # alternative to kpts (A^-1)
    ecutwfc: float = 60.0  # Ry
    ecutrho: float | None = None  # Ry; QE default (4*ecutwfc) when None
    input_data: dict = Field(default_factory=dict)  # nested/flat pw.x namelists


CalculatorConfig = Annotated[
    EmtCalc | TbliteCalc | XtbCalc | MaceCalc | FhiAimsCalc | QeCalc,
    Field(discriminator="type"),
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
    interval: int = 20
    p_translate: float = 0.5
    p_rotate: float = 0.4
    p_conformer: float = 0.1
    max_translate: float = 0.5  # Å
    max_rotate: float = 30.0    # degrees
    refresh_fragments: bool = False
    refresh_scale: float = 1.2
    seed: int | None = None

    @model_validator(mode="after")
    def _moves_nonzero(self) -> MonteCarloSampling:
        if self.p_translate + self.p_rotate + self.p_conformer <= 0:
            raise ValueError("monte_carlo needs at least one move probability > 0")
        return self


SamplingConfig = Annotated[
    MdSampling | RattleSampling | MonteCarloSampling, Field(discriminator="type")
]


# ------------------------------------------------------------------------- selection
class SelectionConfig(TCModel):
    steps: list[str] = Field(default_factory=lambda: ["physicality", "dedup", "diversity"])
    budget: int | None = 50
    min_distance: float = 0.7  # physicality: min interatomic distance (Angstrom)


# -------------------------------------------------------------------------- labeling
class LabelingConfig(TCModel):
    """DFT labeling of the *selected* frames (distinct from the sampling calc).

    Use a DFT calculator (`fhi_aims`/`qe`) here to label the funnel's output with
    energy/forces/stress (+ dipole/polarizability); the cheap `[calculator]` used
    for sampling is configured separately.
    """

    calculator: CalculatorConfig


# --------------------------------------------------------------------------- dataset
class DatasetConfig(TCModel):
    path: str = "dataset"
    format: Literal["extxyz"] = "extxyz"


# --------------------------------------------------------------------- orchestration
class SlurmStage(TCModel):
    """Per-stage Slurm/Apptainer overrides (merged onto built-in defaults)."""

    image: str | None = None  # .sif to `apptainer exec` (default depends on stage)
    partition: str | None = None
    nodes: int = 1
    ntasks: int | None = None
    cpus_per_task: int | None = None
    gpus: int | None = None  # >0 adds --gpus and the apptainer `--nv` flag
    time: str = "01:00:00"
    mem: str | None = None
    qos: str | None = None
    # Override the global runtime/MPI plugin for this stage (None = inherit).
    runtime: Literal["apptainer", "native"] | None = None
    mpi: Literal["pmix", "pmi2", "cray_shasta", "none"] | None = None
    extra_sbatch: list[str] = Field(default_factory=list)  # raw extra "#SBATCH ..." lines
    pre_commands: list[str] = Field(default_factory=list)  # shell lines before the exec
    env: dict[str, str] = Field(default_factory=dict)  # extra exports for this stage


class SlurmConfig(TCModel):
    account: str | None = None
    sif_dir: str = "."  # directory holding the .sif images
    modules: list[str] = Field(default_factory=list)  # `module load` names
    binds: list[str] = Field(default_factory=list)  # apptainer --bind paths
    # How to reach the binaries: our Apptainer images, or binaries already on the
    # host (site modules / conda / EasyBuild). `native` skips the container wrapper.
    runtime: Literal["apptainer", "native"] = "apptainer"
    # Slurm MPI plugin for the multi-node DFT launch. `pmix` for InfiniBand+Slurm
    # clusters (e.g. Leonardo); `cray_shasta` for Cray/Slingshot (e.g. LUMI, which
    # has no pmix); `pmi2` as a portable fallback. Check `srun --mpi=list` on site.
    mpi: Literal["pmix", "pmi2", "cray_shasta", "none"] = "pmix"
    # DFT images used when runtime="apptainer".
    aims_image: str = "traincraft-dft.sif"
    qe_image: str = "traincraft-qe.sif"
    # Full overrides for the DFT engine commands injected into the label stage
    # (container-agnostic plugins). `{binds}`/`{sif_dir}`/`{mpi}` are substituted;
    # when unset the command is composed from runtime + mpi + the images above.
    aims_command: str | None = None
    pw_command: str | None = None
    pre_commands: list[str] = Field(default_factory=list)  # shell lines before every exec
    env: dict[str, str] = Field(default_factory=dict)  # global exports
    stages: dict[str, SlurmStage] = Field(default_factory=dict)


class OrchestrationConfig(TCModel):
    engine: Literal["local", "slurm"] = "local"
    slurm: SlurmConfig | None = None

    @model_validator(mode="after")
    def _slurm_needs_config(self) -> OrchestrationConfig:
        if self.engine == "slurm" and self.slurm is None:
            self.slurm = SlurmConfig()
        return self


# ------------------------------------------------------------------------------ root
class TrainCraftConfig(TCModel):
    run: RunConfig = Field(default_factory=RunConfig)
    geometry: GeometryConfig | None = None
    calculator: CalculatorConfig | None = None
    sampling: SamplingConfig | None = None
    selection: SelectionConfig | None = None
    labeling: LabelingConfig | None = None
    dataset: DatasetConfig | None = None
    orchestration: OrchestrationConfig | None = None
