"""Where the tools of TrainCraft reside

From here, the different modules for specific tasks are called
"""
# Import standard libraries
# import sys
# import os
# import logging
# from pathlib import Path

# Import third-party libraries
# import random
# import numpy as np

# ASE
# from ase import Atoms
# from ase.calculators.espresso import Espresso
# from ase.build import nanotube
# from ase.build import molecule
# from ase.constraints import FixAtoms
# from ase.build.attach import attach
# from ase.io import write
# from ase.io import read
import sys
from ase.visualize import view

# import TrainCraft modules
from config import config
import gengeom

# import calculeitors
import samplers


#   M: gengeom
def generate_geometry():
    mode = config.geom_generation
    system = config.system
    if system == "filled_cnt":
        if mode == "manual":
            cnt = gengeom.create_nanotube()
            molecules = gengeom.add_molecules(cnt)
            system = cnt + molecules
        elif mode == "auto":
            system = gengeom.generate_geometry_filled_cnt_auto()
        gengeom.set_cell(system)
    elif system == "molecular":
        if mode == "manual":
            # TODO
            pass
        elif mode == "auto":
            system = gengeom.generate_geometry_molecular_auto()
    else:
        print(f"ERROR: system {system} not recognised")
        print("""  Options are:
        - cnt
        - filled_cnt
        - molecular
          """)
        sys.exit(1)

    return system


def sample_geometries(system, method="md"):
    if config:
        method = config.sampling_method

    if method == "md":
        samplers.sampling_from_MD(system)
    if method == "rattle":
        samplers.gen_rattled_geometries(system)


def visualize(system):
    """Uses ASE GUI to visualize the system"""
    if config.visualize:
        view(system, repeat=config.vis_repeat)
