`TrainCraft`  automates the process of generating Datasets for training Machine Learning potentials. At the moment it generates randomly filled carbon nanotubes with small molecules, it generates a set of geometries either sampled from a MD (run with DFTTB via tblite) or rattled using HiPhive and calculates the DFT energies and forces of each sampled geometry.

Simple scripts are also provided to collect the final geometries in an .extxyz file or to generate an ASE database.

---

## Required packages

Install the following packages via pip or conda

1. ASE
2. tomlkit

## Optional Packages (depending on your workflow)

To calculate the DFT forces:

3. qe (Quantum Espresso); Although consider using a compiled version in your machine for production calculations. It is way more efficient (but also painful)
5. fhiaims: Requires a license. Accesible for academic use. 


To generate geometries automatically

5. Mdanalysis
6. mdapackmol-fmt (installation via pip)
  It is a forked version of mdapackmol
    `pip install mdapackmol-fmt`
7. Packmol
8. HiPhive (if rattle method is going to be used)


To preoptimize geometries or run MD

9. tblite and tblite-python
11. MACE (requires PyTorch to be installed)
12. torchani

## Installation process example (Linux)

### Installation of conda/miniconda

Download and [install miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install):

    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh

After the installation is finished run: 

    source ~/miniconda3/bin/activate

Initialize conda: 

    conda init --all

Create and activate a new conda environment

    conda create -n traincraft 
    conda activate traincraft

### Installation of the packages for traincraft

Install the required packages, not all are mandatory, but are recomended. The only mandatory ones are ASE and Tomlkit, but then you are missing the main point of using traincraft. 

    conda install -c conda-forge ase tomlkit packmol hiphive mdanalysis tblite tblite-python 

Also need to install mdapackmol-fmt (this will be integrated inside Traincraft in the future, so you will be able to skip this step):

    pip install mdapackmol-fmt

If you want to use MACE, you will need to install pytorch. If you have a CUDA ready device, you can install pytorch with CUDA 12.2 (is a big installation, so it may take a while):

    pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

Once it finishes, you can install MACE as:

    pip install mace-torch


### Installation of DFT codes

If you want to use Quantum Espresso as your ab-initio code, you can also install it via conda, although ideally would be better to have a fresh compilation on your computing machine:

    conda install -c conda-forge qe

Traincraft soon will be able to use [FHI-AIMS](https://fhi-aims.org/) for calculating energy/forces with DFT. There is no conda package, as it requires a license, altough it is free to use for accademic use, it has a voluntary fee to be paid to help running site, and other things.

---

## Minimum input file

This is an input file you can use to run TrainCraft. For this example you will need to have installed the required packages ASE and tomlkit, plus the optional packages Packmol, tblite, Mdanalysis, mdapackmol-fmt, HiPhive and Quantum Espresso. 

```
# Reference input file in TOML format
title = 'Simple TrainCraft input file'

# ----------------------------
# Building geometry parameters
# ----------------------------
[structures]
n_structures = 5  # Total number of initial random geometries to be generated

[molecules]
molecule = 'CO2'
n_molecules = 3
tolerance = 2   # Minimum intermolecular distance for Packmol to distribute the molecules

[cnt]
# n,m nanotube vectors
cnt_n    = 8
cnt_m    = 0
cnt_l    = 2           # Repetition units of a single CNT

[preoptimize]       # The existence of the Table, sets do_preoptimize to True
fmax = 2                 # Max force in atoms

[sampling]          # The existence of the Table, sets do_sampling to True
method = 'rattle'        # 'md'(*) | 'rattle'; How to sample geometries
calculator = 'tblite'    # 'tblite'(*) | 'xtb' | 'ani'; Calculator to use for the sampling method

  [sampling.rattle]    # For more info, look at HiPhive doc: https://hiphive.materialsmodeling.org/advanced_topics/structure_generation.html
    rattle_structures = 10   # how many rattle geometries at each structure

# ---------------------------
# Calculation parameters
# ---------------------------
[calculator]
calculator = 'qe'          # Use Quantum Espresso as calculator to calculate  DFT energies/forces
calculate_forces = true    # If true: DFT forces will be calculated
write_input_file = true    # Writes the input file regardless forces being calculated

  # Set the QE calculator parameters
  [calculator.qe]
    calculation = 'scf'
    # ecutwfc = 30        # Commented plane-wave wave-function cutoff
    ecutwfc = 50          # plane-wave wave-function cutoff
    ecutrho = 600         # density wave-function cutoff
    conv_thr = 1e-8       # DFT self-consistency convergence
    forc_conv_thr = 1e-3  # Force convergence threshold
    etot_conv_thr = 1e-5  # Total energy convergence threshold
    nstep = 100           # Total ionic steps
    pseudo_dir = "/home/user/pseudos/qe/SSSP_1.1.2_PBE_precision/"
    vdw_corr = "xdm"      # use eXchange Dipole Moment method
    occupations = "smearing"   # Add smearing
    smearing = "cold"     # smearing kind
    degauss = 0.02        # smearing amount
    tprnfor = true        # Print forces
    tstress = true        # Print stress tensor

```

## Full Reference input file

`TrainCraft` uses a single toml file as an input. Below is the full reference of tables and keywords.

```
# Reference input file in TOML format
title = 'Input file example for reference'

# Stared options(*) are the default ones

# If [Table] ([preoptimize] | [sample] | [calculator]) is present it will perform the opearation. 
# If not present, the operation will be skipped

# ---------------------------- # Building geometry parameters
# ----------------------------
geom_generation = 'auto'       # 'manual'  : Manual geometry generation
                               # 'auto'(*) : Packmol will be used
                               # false     : Does not generate geometries

[structures]
n_structures = 5  # Total number of initial random geometries to be generated

[molecules]
molecule = 'CO2'
n_molecules = 3
tolerance = 2   # Minimum intermolecular distance for Packmol to distribute the molecules

# compresion_factor = 1.0  # How close molecules will be inside the CNT, when geom_generation = 'manual'

# Axis order in which rotations will be performed, the same axis can be given
# several times, so several rotations will be performed
# rot_axis = ['x', 'y', 'z'](*)
# rot_axis = ['y', 'x', 'z', 'y']
rot_axis = ['x', 'z', 'y']

# Rotations for each molecule (in degrees)
rot_x = 90
rot_y = {'min' = -25.0, 'max' = 25.0}
rot_z = {'min' = 0, 'max' = 360}

[cnt]
# n,m nanotube vectors
cnt_n    = 8
cnt_m    = 0
cnt_l    = 2           # Repetition units of a single CNT
cnt_bond = 1.45        # Bond lenght in the CNT
cnt_gap  = 4           # Distance betwen CNTs (Ang)
# constraints = 'all'  # Fixes all atoms of the Nanotube.  (Lost if using Packmol at the moment)

[preoptimize]       # The existence of the Table, sets do_preoptimize to True
do_preoptimize = true    # Whether to perform preoptimization or not (Takes priority)
fmax = 2                 # Max force in atoms
calculator = 'tblite'    # 'tblite'(*) | 'xtb' | 'ani'; Calculator to use for the optimization method
max_steps = 0.5          # Max distance to move ions

[sampling]          # The existence of the Table, sets do_sampling to True
do_sampling = true       # Whether to perform sampling or not (Takes priority)
method = 'rattle'        # 'md'(*) | 'rattle'; How to sample geometries
calculator = 'tblite'    # 'tblite'(*) | 'xtb' | 'ani'; Calculator to use for the sampling method

  [sampling.md]
    sampling_interval = 20
    temperature = 500
    md_steps = 1000

  [sampling.rattle]    # For more info, look at HiPhive doc: https://hiphive.materialsmodeling.org/advanced_topics/structure_generation.html
    rattle_method = 'mc'     # 'mc'(*) | 'standard'
    rattle_structures = 10   # how many rattle geometries at each structure
    rattle_std = 0.12

    [sampling.rattle.optimize]   # Whether to preoptimize the rattled structures 
      fmax = 10                  # Max force in atoms
      max_steps = 0.5            # Max distance to move ions

# ---------------------------
# Calculation parameters
# ---------------------------
# Quantum Mechanical Calculator
[calculator]
calculator = 'qe'          # Use Quantum Espresso as calculator to calculate  DFT energies/forces
calculate_forces = false   # If true: DFT forces will be calculated
write_input_file = true    # Writes the input file regardless forces being calculated

  # Set the QE calculator parameters
  [calculator.qe]
    calculation = 'scf'
    # ecutwfc = 30        # Commented plane-wave wave-function cutoff
    ecutwfc = 50          # plane-wave wave-function cutoff
    ecutrho = 600         # density wave-function cutoff
    conv_thr = 1e-8       # DFT self-consistency convergence
    forc_conv_thr = 1e-3  # Force convergence threshold
    etot_conv_thr = 1e-5  # Total energy convergence threshold
    nstep = 100           # Total ionic steps
    pseudo_dir = "~/pseudos/qe/SSSP_1.1.2_PBE_precision/"
    vdw_corr = "xdm"
    occupations = "smearing"   # Add smearing
    smearing = "cold"  # smearing kind
    degauss = 0.02     # smearing amount
    tprnfor = true     # Print forces
    tstress = true     # Print stress tensor
    kpts = [1,1,2]     # If not given, Gamma will be used
    # In principle you can keep puting here any QE keyword you need...

    [calculator.qe.pseudos]  # Not mandatory, but recomended
      C = 'C.pbe-n-kjpaw_psl.1.0.0.UPF'
      O = 'O.pbe-n-kjpaw_psl.0.1.UPF'

# ---------------------------
# Visualization parameters
# ---------------------------
[visualize]
# You wanna see the geometries?
visualize = false
# repeat = [1,1,3]

# ---------------------------
# General configuration stuff
# ---------------------------
[config]
# Random things
random_seed = true
# seed = 42  # If random_seed is false, use this seed instead

# Paths and files
folder_prefix = 'dataset'
# logfile = 'logfile.log'
# outdir = 'dir/where/to/generate/the/outputs'

## - Parallelization. If QE < 7  --> set nproc == number of available processor
## - if QE >= 7 do not use nproc, and QE will use all available processors
# nproc = 4   # Number of MPI processes to run qe as: mpiexec -np nproc ...
```
## Citation

A paper about TrainCraft has not yet been published, but if you find TrainCraft useful in your research, you can cite it using the Zenodo DOI: 10.5281/zenodo.8174842
Also, remember to cite the packages you may be using via TrainCraft.

## Acknowledgment

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 847635, under the project UNA Europa, an alliance of universities FOR the emergence of talent and the development of research CAREERs (UNA4CAREER)
 

