import os

import toolkit as tk
from config import config
import logging

# from ase.parallel import MPI

# from gpaw import GPAW
# from ase.calculators.emt import EMT

# TODO:
# v Generate M different geometries
# ~ Add results to the database (Done in separate script so far)
# o Better geometries generation  --> Packmol?
# o Compare similarity between generated structures?
# o Analyse results
# o Add more calculators

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
    - cnt_bond (float), [1.43] CNT C-C bond lenght
    - cnt_gap: (float), [4] Distance betwen neighbour CNTs (Ang)

OUTPUT:
-------
    - Visualizes the system if requested in config file
    - a QE input file for each generated geometry
    - After running the calculation:
        - QE output files

"""

# comm = MPI.COMM_WORLD
# rank = comm.rank

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

    # EMT is just for testing pourposes, EMT breaks CO2 molecules
    # system.calc = EMT()

    calc = tk.set_calculator_parameters()

    # system.set_calculator(calc_QE)  # Deprecated
    system.calc = calc

    # Write input files, should you want to run calculations manually
    system.calc.write_input(system)

    tk.run_MD(system)

    if (config.calculate_e):
        try:
            logging.info('  Calculating energy...')
            system.get_total_energy()
            logging.info('  Energy calculation finished')
        except:
            logging.error(' QEspresso finished with error. Check output')
    if (config.calculate_f):
        try:
            logging.info('  Calculating forces...')
            system.get_forces()
            logging.info('  Forces calculation finished')
        except:
            logging.error(' QEspresso finished with error. Check output')

    os.chdir(config.cwd)
