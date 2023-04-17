import numpy as np
from pathlib import Path
import random
from config import config

from ase.build import nanotube
from ase.build import molecule
from ase.constraints import FixAtoms
from ase.build.attach import attach
# from ase.build.attach import attach_randomly
# from ase.build import make_supercell, find_optimal_cell_shape
# from ase.io import write
# from ase.io import read
# from ase.optimize import BFGS


def set_seed():
    """Duh!"""
    if not config.random_seed:
        random.seed(config.seed)
        np.random.seed(config.seed)


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
      - Randomly rotated and tilted row of molecules along the Z
        axis of the cnt
    """
    molec = config.molec
    n_molecules = config.n_molecules
    compresion_factor = config.compresion_factor
    add_tilt = config.add_tilt
    tilt_factor = config.tilt_factor

    if not config.random_seed:
        random.seed(config.seed)

    # Will arrange the molecules along the Z axis of the CNT
    molecular_distance_z = (
        cnt.get_cell()[2][2])/(n_molecules) * compresion_factor
    # Orient the first molecule
    molecules = molecule(molec)
    molecules.rotate(90, 'x')
    displace = [0, 0, molecular_distance_z * random.random()]
    molecules.translate(displace)
    molecules.rotate(random.randint(0, 360), 'z')
    if add_tilt:
        molecules.rotate(random.randint(0, tilt_factor)*random.choice([-1, 1]),
                         random.choice(['x', 'y']))

    # Add remaining molecules
    mmol = molecule(molec)
    mmol.rotate(90, 'x')
    for m in range(1, n_molecules):
        # mmol.rotate((360/n_molecules * m), 'z')
        mmol.rotate(random.randint(0, 360), 'z')
        if add_tilt:
            mmol.rotate(random.randint(0, tilt_factor)*random.choice([-1, 1]),
                        random.choice(['x', 'y']))

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
    """ Set the calculator parameters"""

    input_params = config.input_params
    pseudos = config.pseudos
    kpts = tuple(config.kpts)

    return input_params, pseudos, kpts

# -----------------------
# Write some files
# -----------------------

# NOTE: the write command does NOT write the calculation parameters
# in the written file :( write( 'system.in', system, format='espresso-in')

# write('system.in', system, format='espresso-in',
#       pseudopotentials={'C': '', 'O': ''})
# traj_file = 'geopt.traj'


# -------------------------
# Path tools
# -------------------------
def set_calculation_folder():
    import uuid
    items = (str(config.prefix),
             str(config.cnt_n),
             str(config.cnt_m),
             str(config.cnt_l),
             str(config.n_molecules) + str(config.molec),
             )
    print(items)
    d = '_'.join(items)
    p = config.outdir / d
    p.mkdir(exist_ok=True)
    calcdir = p / str(uuid.uuid4()).split('-')[0]

    try:
        calcdir.mkdir(exist_ok=False)
    except FileExistsError:
        set_calculation_folder()

    return calcdir
