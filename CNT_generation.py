import os
from pathlib import Path
from ase.calculators.espresso import Espresso
from ase.visualize import view
import toolkit as tk
from config import config

# from gpaw import GPAW
# from ase.calculators.emt import EMT

# TODO:
# - Generate M different geometries
# - Add results to the database
# - Compare similarity between generated structures
# - Add more calculators

"""
Creates CNT randomly filled with small molecules and calculates their energies
and/or forces.

PARAMETERS (read from a config.toml file)
----------
    molec: (str), ['CO2']|'H2O'|'N2'|....
    n_molecules: (int), [4] number of molecules to fill the CNT with
    compresion_factor: (float) controls how close/far the molecules will be
    add_tilt: (bool), whether to add random tilt in either x or y axis
    tilt_factor: (int), [15] max tilt in degrees (can also be )

    cnt_n, cnt_m: (int, int), [8], [0] n and m nanotube vectors
    cnt_l: (int), [2] Repetition units of a single CNT along its Z axis
    cnt_gap: (float), [4] Distance betwen neighbour CNTs (Ang)

OUTPUT:
-------
    Visualizes the system
    After running the calculation:
        QE input file
        QE output files

"""

# Order the chaos
tk.set_seed()

# -----------------------
# Create the geometries
# -----------------------
cnt = tk.create_nanotube()
molecules = tk.add_molecules(cnt)

# Add CNT to molecular system
system = cnt + molecules
tk.set_cell(system, cnt)

# -----------------------
# Visualize the system
# -----------------------
# This is the initial system, not the optimized
view(system, repeat=[2, 2, 2])

# -----------------------
# Run the calculation
# -----------------------

# EMT is just for testing pourposes, EMT breaks CO2 molecules
# system.calc = EMT()

input_params, pseudos, kpts = tk.get_calculator_parameters()
calc_QE = Espresso(input_data=input_params,
                   pseudopotentials=pseudos,
                   kpts=kpts)
system.set_calculator(calc_QE)


os.chdir(tk.set_calculation_folder())
if (config.calculate_e):
    print('Calculating energy...')
    system.get_total_energy()
if (config.calculate_f):
    print('Calculating forces...')
    # system.get_forces()
