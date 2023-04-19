import numpy as np
import sys
from pathlib import Path
import random
from config import config
import logging

from ase.build import nanotube
from ase.build import molecule
from ase.constraints import FixAtoms
from ase.build.attach import attach
from ase.visualize import view

from ase.parallel import parallel_function
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
    logging.info(f'seed : {seed}')


def get_xy_distance(a, b):
    """
    Calculates the distance between 2 vectors in their projections on
    the XY plane (ignores the Z coordinate)
    """
    return np.linalg.norm(a[0:2] - b[0:2])


def create_nanotube():
    """Creates a nanotube

    Parameters
    ----------
    n,m and l are taken from the config file

    Returns
    -------
    a carbon nanotube (ase Atom object)
    """

    # CNT
    cnt_n, cnt_m, cnt_l = config.cnt_n, config.cnt_m, config.cnt_l
    cnt = nanotube(cnt_n, cnt_m, length=cnt_l)

    # mask = [atom.symbol == 'C' for atom in cnt]
    constraints = FixAtoms(mask=[atom.symbol == 'C' for atom in cnt])
    cnt.set_constraint(constraints)

    # ------------------------
    # WORK OUT THE CONSTRAINTS
    # ------------------------

    # Apply constraints based on a radial distance to the z axis of the CNT
    # If the distance between the longitudinal axis Z and the x,y coordintes is
    # bigger than 95% of the radius of the the CNT, them fix the atom

    # NOTE: No need to work out the constraints at the moment, because if
    # I add the CNT with the operator '+' it preserve the constraints.
    # I will keep this here in case I will need it in the future

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
    molec = config.molec
    n_molecules = config.n_molecules
    compresion_factor = config.compresion_factor
    rot_x = config.molec_rot_x
    rot_y = config.molec_rot_y
    rot_z = config.molec_rot_z

    rot_axis = config.molec_rot_axis

    def get_rotation_angle(angle):
        if angle is None:
            return 0
        if isinstance(angle, dict):
            return random.uniform(angle['min'], angle['max'])
        else:
            return angle

    def rotate_molecule(molecule, axis=['x', 'y', 'z']):
        for ax in axis:
            if ax == 'x' and rot_x:
                molecule.rotate(get_rotation_angle(rot_x), 'x')
            if ax == 'y':
                molecule.rotate(get_rotation_angle(rot_y), 'y')
            if ax == 'z':
                molecule.rotate(get_rotation_angle(rot_z), 'z')

    # Will arrange the molecules along the Z axis of the CNT
    molecular_distance_z = (
        cnt.get_cell()[2][2])/(n_molecules) * compresion_factor

    molecules = molecule(molec)

    # Orient the first molecule
    rotate = rot_x or rot_y or rot_z
    if rotate:
        # Rotatmoleculese
        rotate_molecule(molecules, rot_axis)

    # TODO: Add parameters to contrlo displacement
    if True:
        displace = [0, 0, molecular_distance_z * random.random()]
        molecules.translate(displace)

    # Add remaining molecules
    for m in range(1, n_molecules):
        mmol = molecule(molec)
        if rotate:
            rotate_molecule(mmol, rot_axis)
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
    # TODO: If attach_randomly was used check all molecules are inside the CNT
    # (maybe we do not need all molecules inside?)

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
    kpts = tuple(config.kpts)
    if not config.nproc or not isinstance(config.nproc, int) or config.nproc < 1:
        print("Error: invalid value for nproc in config file.")
        sys.exit(1)
    command = f"mpiexec -np {config.nproc} pw.x < espresso.pwi > espresso.pwo"

    return input_params, pseudos, kpts, command

# -----------------------
# Write some files
# -----------------------

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
             str(config.cnt_l),
             str(config.n_molecules) + str(config.molec),
             )
    d = '_'.join(items)
    p = config.outdir / d
    p.mkdir(exist_ok=True)
    calcdir = p / str(uuid.uuid4()).split('-')[0]

    try:
        calcdir.mkdir(exist_ok=False)
    except FileExistsError:
        set_calculation_folder()

    logging.info(f'folder created: {calcdir}')

    return calcdir


def visualize(system):
    if config.visualize:
        view(system, repeat=config.vis_repeat)
        print(config.vis_repeat)
        print(type(config.vis_repeat))
