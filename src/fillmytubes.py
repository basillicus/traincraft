import os

import toolkit as tk
from config import config
import logging

# from ase.parallel import MPI

# from gpaw import GPAW
# from ase.calculators.emt import EMT

# TODO:
# x Generate M different geometries
# - Add results to the database
# - Compare similarity between generated structures
# - Analyse results
# - Add more calculators

"""
Creates CNT randomly filled with small molecules and calculates their energies
and/or forces.

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
    - cnt_gap: (float), [4] Distance betwen neighbour CNTs (Ang)

OUTPUT:
-------
    Visualizes the system if requested in config file
    After running the calculation:
        QE input file
        QE output files

"""

# comm = MPI.COMM_WORLD
# rank = comm.rank

# Order the chaos
tk.set_seed()

structures = config.n_structures
for structure in range(structures):

    # -----------------------
    # Create the geometries
    # -----------------------
    cnt = tk.create_nanotube()
    molecules = tk.add_molecules(cnt)

    # Add molecules and CNT together
    system = cnt + molecules
    tk.set_cell(system, cnt)

    # -----------------------
    # Visualize the system
    # -----------------------
    # This is the initial system, not the optimized
    tk.visualize(system)

    # -----------------------
    # Run the calculation
    # -----------------------

    # EMT is just for testing pourposes, EMT breaks CO2 molecules
    # system.calc = EMT()

    calc = tk.set_calculator_parameters()
    # system.set_calculator(calc_QE)  # Deprecated
    system.calc = calc

    # if rank == 0:
    os.chdir(tk.set_calculation_folder())
    # Write input files, should you want to run calculations manually
    system.calc.write_input(system)

    if (config.calculate_e):
        logging.info('  Calculating energy...')
        system.get_total_energy()
        logging.info('  Energy calculation finished')
    if (config.calculate_f):
        logging.info('  Calculating forces...')
        system.get_forces()
        logging.info('  Forces calculation finished')
    os.chdir(config.cwd)