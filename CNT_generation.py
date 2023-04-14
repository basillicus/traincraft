import numpy as np
import random
from ase.build import nanotube
from ase.build import molecule
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.build.attach import attach, attach_randomly
# from ase.build import make_supercell, find_optimal_cell_shape
# from ase.io import write
# from ase.io import read
# from gpaw import GPAW
# from ase.optimize import BFGS
# from ase.calculators.emt import EMT
from ase.calculators.espresso import Espresso

"""
Creates CNT randomly filled with small molecules and calculates their energies
and/or forces.

PARAMETERS (hardcoded in the script)
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
# ---------------------------
# System Geometry parameters
# ---------------------------

# TODO: Read the configuration from a config file?
# including the input_params (instead of hardcoding it in the script)

# Molecules
molec = 'CO2'
n_molecules = 4
compresion_factor = 1  # How close the molecules will be inside the CNT
add_tilt = True
tilt_factor = 15    # max randomly molecular tilt around either x or y

# CNT
cnt_n, cnt_m = 8, 0  # n,m nanotube vectors
cnt_l = 2    # Repetition units of a single CNT
cnt_gap = 4  # Distance betwen CNTs (Ang)

# Random things
random_seed = True
if not random_seed:
    seed = 42
    random.seed(seed)

# ---------------------------
# Calculation parameters
# ---------------------------

input_params = {
    "ecutwfc": 45,     # plane-wave wave-function cutoff
    "ecutrho": 180,    # density wave-function cutoff,
    "conv_thr": 1e-6,  # DFT self-consistency convergence
    "pseudo_dir": "/home/david/pseudos/qe/SSSP_1.1.2_PBE_precision/",
    "vdw_corr": "xdm",
    "occupations": 'smearing',   # Add smearing
    "smearing": 'cold',  # smearing kind
    "degauss": 0.02,     # smearing amount
    "tprnfor": True,     # Print forces
    "tstress": True      # Print stress tensor
}
# Add a new paramter to the input_params dict
input_params["calculation"] = 'scf'

# define the pseudopotentials
pseudos = {"C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
           "O": "O.pbe-n-kjpaw_psl.0.1.UPF"}

kpts = (1, 1, 1)
# kpts = (2, 2, 2)

# -----------------
# Build the System
# -----------------


def get_xy_distance(a, b):
    """
    Calculates the distance between 2 vectors in their projections on
    the XY plane (ignores the Z coordinate)
    """
    return np.linalg.norm(a[0:2] - b[0:2])


def create_nanotube():
    # CNT
    cnt = nanotube(cnt_n, cnt_m, length=cnt_l)

    # mask = [atom.symbol == 'C' for atom in cnt]
    # constraints = FixAtoms(mask=[atom.symbol == 'C' for atom in cnt])
    constraints = FixAtoms(mask=[atom.symbol == 'C' for atom in cnt])
    cnt.set_constraint(constraints)

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

    # Will arrange the molecules along the Z axis of the CNT
    molecular_distance_z = (
        cnt.get_cell()[2][2])/(n_molecules) * compresion_factor
    print(molecular_distance_z)
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
    """Set the cell parameter of the system according to the CNT dimensions"""
    # TODO: If attach_randomly was used check all molecules are inside the CNT
    # (maybe we do not need all molecules inside?)

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


cnt = create_nanotube()
molecules = add_molecules(cnt)

# Add the CNT to the molecular system
# system = attach(molecules,cnt, distance=1.9 )  # Does not preserve constraints
system = cnt + molecules      # It does preserve constraints
set_cell(system, cnt)


# ------------------------
# WORK OUT THE CONSTRAINTS
# ------------------------

# Apply constraints based on a radial distnnce to the z axis of the CNT
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


# -----------------------
# Write some files
# -----------------------

# NOTE: the write command does NOT write the calculation parameters
# in the written file :( write( 'system.in', system, format='espresso-in')

# write('system.in', system, format='espresso-in',
#       pseudopotentials={'C': '', 'O': ''})
# traj_file = 'geopt.traj'

# -------------------------
# Visualize the system
# -------------------------

# This is the initial system, not the optimized
view(system, repeat=[2, 2, 2])

# ------------------------
# Run the calculation
# ------------------------

# TODO: Run it with GPAW, FHI-AIMS...?
#
# EMT is just for testing pourposes, in fact the CO2 molecules break
# system.calc = EMT()
calc_QE = Espresso(input_data=input_params,
                   pseudopotentials=pseudos,
                   kpts=kpts)
system.set_calculator(calc_QE)

# print('Calculating energy...')
# system.get_total_energy()
# print('Calculating forces...')
# system.get_forces()

# TODO: Add results to the database
# Generate M different geometries
# Compare similarity between generated molecules
