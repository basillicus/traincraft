import sys
import logging
# from pathlib import Path

import random
import numpy as np

from config import config

from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.build import nanotube
from ase.build import molecule
from ase.constraints import FixAtoms
from ase.build.attach import attach
from ase.visualize import view

# from ase.parallel import parallel_function
# from ase.build.attach import attach_randomly
# from ase.build import make_supercell, find_optimal_cell_shape
# from ase.io import write
# from ase.io import read
# from ase.optimize import BFGS


def set_seed():
    """Duh!"""
    if config.random_seed:
        seed = random.randrange(2**32 - 2)
    else:
        seed = config.seed
    random.seed(seed)
    np.random.seed(seed)
    logging.info(f'seed: {seed}')


def get_xy_distance(a, b):
    """
    Calculates the distance between 2 vectors in their projections on
    the XY plane (ignores the Z coordinate)
    """
    return np.linalg.norm(a[0:2] - b[0:2])


def create_nanotube():
    """Clue: Creates a ...

    Parameters
    ----------
    n,m and lenght are taken from the config file

    Returns
    -------
    a carbon nanotube (ase Atom object). Do not eat!
    """

    n, m, lenght = config.cnt_n, config.cnt_m, config.cnt_l
    bond_lenght = config.cnt_bond
    constraints = config.cnt_constraints

    cnt = nanotube(n, m, length=lenght, bond=bond_lenght)

    # mask = [atom.symbol == 'C' for atom in cnt]
    if constraints == 'all':
        constraints = FixAtoms(mask=[atom.symbol == 'C' for atom in cnt])
        cnt.set_constraint(constraints)

    # ------------------------
    # WORK OUT THE CONSTRAINTS
    # ------------------------
    # NOTE: Sorry for the mess below... A more flexible way of adding
    # constraints will be useful, but it may take some time before is
    # implemented. I keep this here for the future...

    # Apply constraints based on a radial distance to the z axis of the CNT
    # If the distance between the longitudinal axis Z and the x,y coordintes is
    # bigger than 95% of the radius of the the CNT, then fix the atom

    # system_cto = cnt.get_center_of_mass()       # center of Mass of the CNT
    # cnt_atom =  cnt.get_positions()[1]          # take a random atom of the CNT
    # r = get_xy_distance(system_cto, cnt_atom)   # radi of the CNT

    #   loop over atom:        _/---------------------------------------\_
    # mask = [get_xy_distance(atom.position, system_cto) >= r*0.95 for atom in system]
    #         |------------------condicion-----------------------|
    #
    # constraints = FixAtoms(mask=mask)
    # system.set_constraint(constraints)

    return cnt


def add_molecules(cnt):
    """ Add ONLY the molecules along the Z axis of the nanotube
    Does not add the nanotube itself

    Parameters
    ----------
      - cnt: nanotube

    Returns
    -------
      -  Molecules along the Z axis of the cnt, rotated and displaced according
         to the input parameters
    """
    # Get the parameters from config object
    molec = config.molec
    n_molecules = config.n_molecules
    compresion_factor = config.compresion_factor
    rot_x = config.molec_rot_x
    rot_y = config.molec_rot_y
    rot_z = config.molec_rot_z

    rot_axis = config.molec_rot_axis

    def rotate_molecule(molecule, axis=['x', 'y', 'z']):
        for ax in axis:
            if ax == 'x':
                molecule.rotate(get_rotation_angle(rot_x), 'x')
            if ax == 'y':
                molecule.rotate(get_rotation_angle(rot_y), 'y')
            if ax == 'z':
                molecule.rotate(get_rotation_angle(rot_z), 'z')

    def get_rotation_angle(angle):
        if angle is None:
            return 0
        if isinstance(angle, dict):
            return random.uniform(angle['min'], angle['max'])
        else:
            return angle

    # It will arrange the molecules along the Z axis of the CNT
    if n_molecules > 0:
        molecular_distance_z = (
            cnt.get_cell()[2][2])/(n_molecules) * compresion_factor

    rotate = rot_x or rot_y or rot_z

    molecules = Atoms()
    # Add molecules
    for m in range(0, n_molecules):
        mmol = molecule(molec)
        if rotate:
            rotate_molecule(mmol, rot_axis)

        if m == 0:
            # TODO: Add parameters to control displacement
            if True:
                # Displace the first molecule
                displace = [0, 0, molecular_distance_z * random.random()]
                mmol.translate(displace)
            molecules = mmol
        else:
            molecules = attach(molecules, mmol,
                               distance=molecular_distance_z,
                               direction=(0, 0, 1))
    return molecules


def set_cell(system, cnt):
    """Set the cell parameter of the system according to the CNT dimensions

    Parameters
    ----------
    System: The CNT with the molecules
    cnt: The CNT itself
    """
    # TODO: If attach_randomly was used check that all molecules are
    # inside the CNT (maybe we do not need all molecules inside?)

    cnt_gap = config.cnt_gap

    # Get the limits of the CNT
    cnt_xyz = cnt.get_positions()
    min_x, max_x = min(cnt_xyz[:, 0]), max(cnt_xyz[:, 0])
    min_y, max_y = min(cnt_xyz[:, 1]), max(cnt_xyz[:, 1])
    # min_z, max_z = min(cnt_xyz[:, 2]), max(cnt_xyz[:, 2])
    cell_z = cnt.get_cell()[2][2]

    cell_x = [(max_x - min_x + cnt_gap), 0, 0]
    cell_y = [0, (max_y - min_y + cnt_gap), 0]
    cell_z = cnt.get_cell()[2]  # Use the lattice parameter Z of the CNT
    system.set_cell([cell_x, cell_y, cell_z])


def get_calculator_parameters():
    """ Get the calculator parameters from config object"""

    input_params = config.input_params
    pseudos = config.pseudos
    if config.kpts is not None:
        kpts = tuple(config.kpts)
    else:
        kpts = None
    nproc = config.nproc

    if nproc is not None and (not isinstance(nproc, int) or nproc < 1):
        print("Error: invalid value for nproc in config file.")
        sys.exit(1)

    command = None
    if nproc is not None:
        command = f"mpiexec -np {nproc} pw.x < espresso.pwi > espresso.pwo"
    return input_params, pseudos, kpts, command


def set_calculator_parameters():
    """ Set the calculator parameters from config object"""

    # TODO: When adding support for more calculators, come here and tweak it
    input_params, pseudos, kpts, command = get_calculator_parameters()
    if command:
        calc = Espresso(input_data=input_params,
                        pseudopotentials=pseudos,
                        kpts=kpts,
                        command=command)

    else:
        calc = Espresso(input_data=input_params,
                        pseudopotentials=pseudos,
                        kpts=kpts)
    return calc
# -------------------------
# Write files
# -------------------------

# TODO: write input files without running the calculation?
# Is it possible with ASE?

# NOTE: the write command does NOT write the calculation parameters
# in the written file :( write( 'system.in', system, format='espresso-in')

# write('system.in', system, format='espresso-in',
#       pseudopotentials={'C': '', 'O': ''})
# traj_file = 'geopt.t eraj'


# -------------------------
# Path tools
# -------------------------
# @ parallel_function
# def set_calculation_folder(parallel=False):
def set_calculation_folder():
    import uuid
    items = (str(config.prefix),
             str(config.cnt_n) + '-' + str(config.cnt_m),
             'len-' + str(config.cnt_l),
             str(config.n_molecules) + str(config.molec),
             )
    d = '_'.join(items)
    p = config.outdir / d
    p.mkdir(exist_ok=True)
    calcdir = p / str(uuid.uuid4()).split('-')[0]

    try:
        calcdir.mkdir(exist_ok=False)
    except FileExistsError:
        # Really??! o.O Well, try again...
        set_calculation_folder()

    logging.info(f'folder created: {calcdir}')

    return calcdir


def visualize(system):
    """Uses ASE GUI to visualize the system"""
    if config.visualize:
        view(system, repeat=config.vis_repeat)
