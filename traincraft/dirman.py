import logging
import os

import calculeitors
from config import config

"""DIRectory MANager: Module to handle paths, folders, files and all that jazz"""

# TODO: Implement and use this module
SAMPLING_FOLDER = 'sampled_geometries'
DFT_FOLDER = 'dft_forces'

def write_input_file(system):
    # Copy the original system to not modify it when assigning the calculator
    system_copy = system
    # This calculator is attached  only to save an input file
    calc = calculeitors.set_DFT_calculator_parameters()
    system_copy.calc = calc
    # Write input files, should you want to run calculations manually
    system_copy.calc.write_input(system_copy)
    # KKK
    system.calc.write_input(system)

# -------------------------
# Path tools
# -------------------------
def set_calculation_folder():
    import uuid
    # TODO: Change items according to if it is a crystal, surface, cnt...
    items = (str(config.prefix),
             str(config.cnt_n) + '-' + str(config.cnt_m),
             'len-' + str(config.cnt_l),
             str(config.n_molecules) + str(config.molec),
             )
    d = '_'.join(items)
    p = config.outdir / d
    os.makedirs(p, exist_ok=True)

    calcdir = p / str(uuid.uuid4()).split('-')[4]

    try:
        calcdir.mkdir(exist_ok=False)
    except FileExistsError:
        # Really??! o.O Well, try again...
        set_calculation_folder()
    logging.info(f'folder created: {calcdir}')

    return calcdir
