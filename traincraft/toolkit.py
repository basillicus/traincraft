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
from ase.visualize import view

# import TrainCraft modules
from config import config
import gengeom
import calculeitors


#   M: gengeom
def generate_geometry():

    mode = config.geom_generation
    if mode == 'manual':
        cnt = gengeom.create_nanotube()
        molecules = gengeom.add_molecules(cnt)
        system = cnt + molecules
    elif mode == 'auto':
        system = gengeom.generate_geometry_with_packmol()

    gengeom.set_cell(system)

    return system


def sample_geometries(system, method='md'):
    if config:
        method = config.sampling_method

    if method == 'md':
        calculeitors.sampling_from_MD(system)
    if method == 'rattle':
        gengeom.gen_rattled_geometries(system)


def visualize(system):
    """Uses ASE GUI to visualize the system"""
    if config.visualize:
        view(system, repeat=config.vis_repeat)
