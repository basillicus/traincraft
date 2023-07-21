import os
import logging

from config import config
import toolkit as tk
import calculeitors
import utils
import dirman
# import gengeom

# TODO:
# v Generate M different geometries
# ~ Add results to the database (Done in separate script so far) v Better geometries generation  --> Packmol
# v Generate sampling geometries --> From Tight Bidning MD
# o ALLOW FOR FLEXIBLE WORKFLOW:
#   - Sometimes you want generate geometries
#   - sometimes you want to go through the generated geometries and do sth
# o Add gneration of different systems: Molecules, molecs on surfaces, crystals, intercalated molecules...
# o Analyse results
# o Add more calculators
# o Compare similarity between generated structures? Not at the moment...

"""
Creates CNT randomly filled with small molecules and calculates their
DFT energies and forces, generating a dataset for training ML potentials.

Parameters
----------
(read from a config.toml file)
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

Output
------
    - Visualizes the system if requested in config file
    - a QE input file for each generated geometry
    - md.traj trajectories
    - sampled geometries in .extxyz format
    - After running the QM calculation:
        - QE output files

"""
# TODO: Set a flexible (not lineal) WORKFLOW
# # Workflow = ['preoptparams', 'gen_geom', 'preopt', 'sample_geometries', ...]
# for work in config.workflow:
#   call work(params)
#
def run_traincraft():
    """Runs the main logic of code"""
    structures = config.n_structures
    for structure in range(structures):
        # Go to the required folder
        os.chdir(dirman.set_calculation_folder())

        # -----------------------
        # Create the geometries
        # -----------------------

        system = None
        if config.geom_generation:
            system = tk.generate_geometry()

        # -----------------------
        # Visualize the system
        # -----------------------
        # This is the initial system, not the optimized
        if config.visualize and system is not None:
            tk.visualize(system)

        # -----------------------
        # Run the preoptimization
        # -----------------------
        if config.do_preoptimize and system is not None:
            calculeitors.preotimize(system)

        # -----------------------
        # Sample geometries
        # -----------------------
        # Generate rattled structures via MD
        if config.do_sampling and system is not None:
            tk.sample_geometries(system)

        # -----------------------
        # Run the QM calculations
        # -----------------------
        if config.calculator is not None:
            # If the calculator has been set, write the input file
            # Save initial Input Files
            if config.calculator_writeInput and system is not None:
                dirman.write_input_file(system)
            if config.calculate_f:
                if config.do_sampling:
                    # Go through all generated .extxyz files and calcute their DFT forces
                    if config.sampling_method == 'md':
                        calculeitors.get_QM_forces_from_sampled(sampled_by='md')
                    if config.sampling_method == 'rattle':
                        calculeitors.get_QM_forces_from_sampled(sampled_by='rattle')
                else:
                    # If not sampled geometries, just get the forces of the generated geometry
                    calculeitors.get_QM_forces(system)

        # Go back and repeat
        os.chdir(config.cwd)
