`TrainCraft`  automatizes the process of generating Datasets for training Machine Learning potentials. At the moment it generates randomly filled carbon nanotubes with small molecules, it generates a set of geometries either sampled from a MD (run with DFTTB via tblite) or rattled using HiPhive and calculates the DFT energies and forces of each sampled geometry.

Few simples scripts are also provided, to collect all the final geometries in an .extxyz file or to generate an ASE database.


TODO: 
 - At the moment it only generates filled nanotubes, but it should be easily extensible to other systems
 - Generate Machine Larned force fields for the required system 
 - Run molecular dynamics with the generated force fields
 - Calculate the Raman spectra from the MD using a ML methodology

---

## Required packages

Install the following packages via pip or conda

1. ASE
2. tomlkit

## Optional Packages (depending on your workflow)

To calculate the DFT forces:
3. qe (Quantum Espresso); Although consider using a compiled version in your machine for production calculations. It is way more efficient (but also painful)

To generate geometries automatically
4. Mdanalysis
5. mdapackmol-fmt (installation via pip)
  It is a forked version of mdapackmol
    `pip install mdapackmol-fmt`
6. Packmol
7. HiPhive (if rattle method is going to be used)


To preoptimize geometries or run MD
8.tblite and tblite-python
9.torchani

---

## Input file

`Fillmytubes` uses a toml file as input. As an example, use:

```
## Keywords

title = 'Input file example. Use it as a reference'

# ------------------------------
# Building geometry parameters
# ------------------------------

# By default, geometries will be generated using Packmol code. 
# For manual control, set it to manual. Then rot_axis and rot_* will be used.
# geom_generation = 'manual'    # Default: 'packmol'


[structures]
n_structures = 10  # Total number of random geometries to be generated

[molecules]
molecule = 'CO2'
n_molecules = 6
tolerance = 2.0  # Minimum distance between molecules when packing (for Packmol)

# If geom_generation == 'manual':

    compresion_factor = 1.0  # How close molecules will be inside the CNT

    # Axis order in which molecule rotations will be performed, the same axis can be given    # several times, so several rotations will be performed.
    # Default is ['x', 'y', 'z']
    # rot_axis = ['y', 'x', 'z', 'y']
    rot_axis = ['x', 'z', 'y']
    
    # Rotations for each molecule (in degrees) 
    rot_x = 90
    rot_y = {'min' = -25.0, 'max' = 25.0}
    rot_z = {'min' = 0, 'max' = 360}

[cnt]
# n,m nanotube vectors
cnt_n = 8
cnt_m = 0
cnt_l = 2        # Repetition units of a single CNT
cnt_bond = 1.42  # Bond lenght in the CNT
cnt_gap = 4      # Distance betwen CNTs (Ang)
# constraints = 'all'  # all: fix all nantoube atoms | none: none are fix  

# ---------------------------
# Calculation parameters
# ---------------------------
[calculator]
calculator = 'qe'         # Use Quantum Espresso to calculate enrgies/forces
calculate_forces = false  # If true: DFT forces are calculated from MD sampled geometries

  # Set the calculator parameters
  [calculator.qe]
    calculation = 'scf'    # QE kind of calculation: scf|relax|vc-relax
    # calculation = 'relax'  # Commented calculation    
    # ecutwfc = 30         # commented plane-wave wave-function cutoff
    ecutwfc = 45           # plane-wave wave-function cutoff
    ecutrho = 180          # density wave-function cutoff,
    conv_thr = 1e-8        # DFT self-consistency convergence
    forc_conv_thr = 1e-3   # Force convergence threshold
    etot_conv_thr = 1e-5   # Total energy convergence threshold
    nstep = 100            # Total ionic steps
    pseudo_dir = "/home/user/pseudos/qe/SSSP_1.1.2_PBE_precision/"
    vdw_corr = "xdm"       # Dispersion interactions
    occupations = "smearing"   # Add smearing
    smearing = "cold"      # smearing kind
    degauss = 0.02         # smearing amount
    tprnfor = true         # Print forces
    tstress = true         # Print stress tensor
    kpts = [1,1,2]         # [Kx, Ky, Kz] or Not given|None: will use Gamma point only

    [calculator.qe.pseudos]
      C = 'C.pbe-n-kjpaw_psl.1.0.0.UPF'
      O = 'O.pbe-n-kjpaw_psl.0.1.UPF'

# ---------------------------
# Visualization parameters
# ---------------------------
[visualize]
# You wanna see the geometries?
visualize = true
# repeat = [1,1,3]   # Repetion of the unit cell for visualization pourposes only

# ---------------------------
# General configuration stuff
# ---------------------------
[config]
# Random things
random_seed = true
# seed = 42  # If random_seed is false, use this seed instead

# Paths and files
folder_prefix = 'dataset'
logfile = 'logfile.log'
# outdir = './'

# Parallelization (Experimental!)
# nproc = 4   # Number of MPI processes to run qe as: mpiexec -np nproc ...
```
