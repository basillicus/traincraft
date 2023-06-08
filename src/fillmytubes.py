import os

import toolkit as tk
from config import config

# import logging

# from ase.parallel import MPI

# from gpaw import GPAW
# from ase.calculators.emt import EMT

# TODO:
# v Generate M different geometries
# ~ Add results to the database (Done in separate script so far) v Better geometries generation  --> Packmol
# v Generate sampling geometries --> From Tight Bidning MD
# o ALLOW FOR FLEXIBLE WORKFLOW:
#   - Sometimes you want generate geometries
#   - sometimes you want to go through the generated geometries and do sth
# o Analyse results
# o Add more calculators
# o Compare similarity between generated structures? Not at the moment...

"""
Creates CNT randomly filled with small molecules and calculates their DFT energies
and forces, generating a dataset for training ML potentials.

PARAMETERS (read from a config.toml file)
----------
    - molec: (str), ['CO2']|'H2O'|'N2'|....
    - n_molecules: (int), [4] number of molecules to fill the CNT with
    - compresion_factor: (float) controls how close/far the molecules will be
    - rot_?: rotation along the ? axis (x, y or z). Can be either a fixed value
      or a dict with min and max values. If not given, will be None
        rot_x = 78.3   or rot_y = {'min' = -12 , 'max' = 32.3}
    - rot_axis: List with the order of axes rotations will be performed on:
        rot_axis = ['x', 'y', 'z'] (default)
        rot_axis = ['y', 'z', 'x', 'y']

    - cnt_n, cnt_m: (int, int), [8], [0] n and m nanotube vectors
    - cnt_l: (int), [2] Repetition units of a single CNT along its Z axis
    - cnt_bond (float), [1.43] CNT C-C bond lenght
    - cnt_gap: (float), [4] Distance betwen neighbour CNTs (Ang)

OUTPUT:
-------
    - Visualizes the system if requested in config file
    - a QE input file for each generated geometry
    - md.traj trajectries
    - sampled geometries in .extxyz format
    - After running the QM calculation:
        - QE output files

"""
#
# Order the chaos
tk.set_seed()

structures = config.n_structures
for structure in range(structures):
    # if rank == 0:
    os.chdir(tk.set_calculation_folder())

    # -----------------------
    # Create the geometries
    # -----------------------

    system = tk.generate_geometry()

    # Add molecules and CNT together
    tk.set_cell(system)

    # -----------------------
    # Visualize the system
    # -----------------------
    # This is the initial system, not the optimized
    tk.visualize(system)

    # -----------------------
    # Run the calculation
    # -----------------------

    calc = tk.set_DFT_calculator_parameters()
    # system.set_calculator(calc_QE)  # Deprecated
    system.calc = calc
    # Write input files, should you want to run calculations manually
    system.calc.write_input(system)

    # Generate rattled structures via MD using DFTB
    # TODO: Allow somehow control this parameteres too in input file
    # But it needs to be a good trade between usability and flexibilty
    tk.run_MD(system, md_steps=5000, fmax=0.05)  # Fmax is used for preoptimization

    # Go through all the generated extxyz files and calcute their DFT forces
    tk.get_QM_forces()

    os.chdir(config.cwd)
