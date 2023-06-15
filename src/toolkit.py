import sys
import os
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

from ase.io import write
from ase.io import read

# Module:
# ------
#   M: utils
def set_seed():
    """Duh!"""
    if config.random_seed:
        seed = random.randrange(2**32 - 2)
    else:
        seed = config.seed
    random.seed(seed)
    np.random.seed(seed)
    logging.info(f'seed: {seed}')


#   M: utils
def get_xy_distance(a, b):
    """
    Calculates the distance between 2 vectors in their projections on
    the XY plane (ignores the Z coordinate)
    """
    return np.linalg.norm(a[0:2] - b[0:2])

# -------------------------
# Write files
# -------------------------

# write('system.in', system, format='espresso-in',
#       pseudopotentials={'C': '', 'O': ''})
# traj_file = 'geopt.traj'

# File manipulation stuff
#   M: utils
def get_parameter_from_pdb(pdb_file, pattern=None):
    with open(pdb_file, 'r') as file:
        for line in file:
            if pattern in line:
                return line

#   M: utils
def insert_line_in_pdb(pdb_file, new_line, insert_index):
    # Read the contents of the PDB file
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    # Insert the new line at the desired index
    lines.insert(insert_index, new_line + '\n')

    # Write the modified lines back to the PDB file
    with open(pdb_file, 'w') as file:
        file.writelines(lines)

# -------------------------
# Path tools
# -------------------------
#   M: utils
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
    p.mkdir(exist_ok=True)
    calcdir = p / str(uuid.uuid4()).split('-')[4]

    try:
        calcdir.mkdir(exist_ok=False)
    except FileExistsError:
        # Really??! o.O Well, try again...
        set_calculation_folder()

    logging.info(f'folder created: {calcdir}')

    return calcdir


# Module:
# ------
#   GenGeom
#   M: gengeom
def generate_geometry():

    mode = config.geom_generation
    if mode == 'manual':
        cnt = create_nanotube()
        molecules = add_molecules(cnt)
        system = cnt + molecules
    elif mode == 'auto':
        system = generate_geometry_with_packmol()

    return system


#   M: gengeom
def generate_geometry_with_packmol():
    """
    Generate geometry using PACKMOL via mdapackmol and MDAnalysis to
    preserve the topology after Packmol
    """
    import MDAnalysis as mda
    import mdapackmol

    # Read required config info
    molec = config.molec
    n_molecules = config.n_molecules
    tolerance = config.tolerance
    cnt_gap = config.cnt_gap
    random_seed = config.random_seed
    # seed = config.seed

    cnt_file = 'cnt_byFillMyTubes.pdb'
    mol_file = 'mol_byFillMyTubes.pdb'

    # Generate the Nanotube and the molecules
    cnt = create_nanotube()
    cnt.center(vacuum=cnt_gap/2, axis=(0, 1))
    write(cnt_file, cnt)
    # write('cnt_byFillMyTubes.xyz', cnt)
    mmol = molecule(molec)
    mmol.center()
    write(mol_file, mmol)
    # write('mol_byFillMyTubes.xyz', mmol)
    cell_parameters = get_parameter_from_pdb(cnt_file, 'CRYST1')

    # Get ranges for PACKMOL
    cnt_xyz = cnt.get_positions()
    min_x, max_x = min(cnt_xyz[:, 0]), max(cnt_xyz[:, 0])
    min_y, max_y = min(cnt_xyz[:, 1]), max(cnt_xyz[:, 1])
    min_z, max_z = min(cnt_xyz[:, 2]), max(cnt_xyz[:, 2])

    cnt_z_length = max_z - min_z
    cnt_x_diameter = (max_x - min_x)/2
    cnt_y_diameter = (max_y - min_y)/2

    cnt_diameter = (cnt_x_diameter + cnt_y_diameter)/2

    # cell_z = cnt.get_cell()[2][2]

    # load individual molecule files into MDAnalysis Universe
    _cnt = mda.Universe(cnt_file)
    # _cnt = mda.Universe('cnt_byFillMyTubes.xyz')

    _molecule = mda.Universe(mol_file)
    # _molecule = mda.Universe('mol_byFillMyTubes.xyz')

    # This is confusing, but if we want to recreate the Packmol geometries we need to pass a
    # different seed each time generated with the seed we gave to FillMyTubes
    if random_seed == False:
        seed = random.randint(1, 2**16-2)
    else:
        seed = -1

    # call Packmol with MDAnalysis objects as arguments
    # the 'instructions' allow for any valid Packmol commands, one per line
    system = mdapackmol.packmol(
        [mdapackmol.PackmolStructure(
            _cnt, number=1,
            instructions=[f'fixed 0. 0. {cnt_z_length/2} 0. 0. 0.', 'center']),
            mdapackmol.PackmolStructure(
            _molecule, number=n_molecules,
            instructions=[f'inside cylinder 0. 0. 0. 0. 0. 1. {cnt_diameter-0.05} {cnt_z_length-tolerance/3}'])],
        tolerance=tolerance,
        seed=seed
    )

    # We do not need MDAnalysis Universe at the moment, we want ASE Atoms
    # so we read the geometry generated by packmol again with ASE

    # First we need to add the cell parameters that are missing after packmol
    insert_line_in_pdb('output.pdb', cell_parameters[:-1], 5)
    # Read the output generated by Packmol
    system = read('output.pdb')

    # Clean up and sort out the mess
    for f in [cnt_file, mol_file, 'packmol.stdout', 'output.pdb']:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

    return system


#   M: gengeom
def create_nanotube():
    """Clue: Creates a ...

    Parameters
    ----------
    n,m and lenght are taken from the config file

    Returns
    -------
    a carbon nanotube (ase Atom object)
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


#   M: gengeom
# TODO:
# Move manipulation functions out of this function so I can use it to manipulate
# individual geometries so we can pass to PACKMOL molecules in the desired orientation.
#   - My be useful when placing molecules outside the nanotube or on surfaces.
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
        print('molecular distance', molecular_distance_z)

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


#   M: gengeom
def set_cell(system):
    """Set the cell parameter of the system

    Parameters
    ----------
    System: (ASE Atoms)

    Returns
    -------
        Returns nothing. It modifies the system itself
    """

    cnt_gap = config.cnt_gap
    system.center(vacuum=cnt_gap/2, axis=(0, 1))

    # OLD Manual Way
    # # Get the limits of the CNT
    # cnt_xyz = cnt.get_positions()
    # min_x, max_x = min(cnt_xyz[:, 0]), max(cnt_xyz[:, 0])
    # min_y, max_y = min(cnt_xyz[:, 1]), max(cnt_xyz[:, 1])
    # # min_z, max_z = min(cnt_xyz[:, 2]), max(cnt_xyz[:, 2])
    # cell_z = cnt.get_cell()[2][2]

    # cell_x = [(max_x - min_x + cnt_gap), 0, 0]
    # cell_y = [0, (max_y - min_y + cnt_gap), 0]
    # cell_z = cnt.get_cell()[2]  # Use the lattice parameter Z of the CNT
    # system.set_cell([cell_x, cell_y, cell_z])


#   M: gengeom
def gen_rattled_geometries(system, min_distance=1.3):
    """Generate rattled geometries using HiPhive

    Parameters
    ----------
    System: (ASE Atoms)

    Generates
    ---------
    Rattled geoemetries: .extxyz files
    """
    # """
    # Generate displaced structures using
    # * Standard rattle, gaussian displacements
    # * Monte Carlo rattle, gaussian displacements w penatly for short atomic dists
    # * Phonon rattle, construct a rough guess of fc2 and populate phonon modes with
    #   thermal energy corresponding to a temperature
    #
    # The hyper parameters for the different methods are chosen such that the
    # magnitude of the displacements will be roughly the same
    #
    # This script may take a few minutes to run.
    # """

    from ase.io import write
    # from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
    # from trainstation import Optimizer
    # from hiphive.utilities import prepare_structures
    # from hiphive.structure_generation import generate_rattled_structures
    from hiphive.structure_generation import (generate_rattled_structures,
                                              generate_mc_rattled_structures,
                                              generate_phonon_rattled_structures)

    method = config.sampling_calculator
    rattled_structures = config.rattle_structures
    rattle_std = config.rattle_std
    rattle_method = config.rattle_method

    calc = get_aproximate_calculator(method)

    # supercell = prim.repeat(size)
    reference_positions = system.get_positions()

    # write('reference_structure.xyz', supercell)

    # standard rattle
    if rattle_method == 'standard':
        logging.info(f'Rattling: {rattle_method}; Rattled structures: {rattled_structures}; Rattle std: {rattle_std}')
        for i in range(rattled_structures):
            seed = np.random.randint(2**32-2)
            structure_rattle = generate_rattled_structures(system, 1, rattle_std, seed=seed)
            # TODO: Work out the folder storage/retrieval
            write('structures_rattle_' + str(i) + '.extxyz', structure_rattle)

    # Monte Carlo rattle
    if rattle_method == 'mc':
        logging.info(f'Rattling: {rattle_method}; Rattled structures: {rattled_structures}; Rattle std: {0.25 * rattle_std}; Min distance: {min_distance}')
        for i in range(rattled_structures):
            seed = np.random.randint(2**32-2)
            try:
                structures_mc_rattle = generate_mc_rattled_structures(
                    system, 1, 0.25*rattle_std, min_distance, n_iter=20, seed=seed)
                # TODO: Work out the folder storage/retrieval
                write(f'structures_mc_rattle_{i}.extxyz', structures_mc_rattle)
            except Exception as e:
                logging.warn(e)
                pass

    #   # Phonon rattle

    #   # initial model
    #   cs = ClusterSpace(supercell, [5.0])
    #   sc = StructureContainer(cs)
    #   rattled_structures = generate_rattled_structures(supercell, 1, 0.05)
    #   rattled_structures = prepare_structures(rattled_structures, supercell, calc)

    #   for structure in rattled_structures:
    #       sc.add_structure(structure)
    #   opt = Optimizer(sc.get_fit_data(), train_size=1.0)
    #   opt.train()
    #   fcp = ForceConstantPotential(cs, opt.parameters)

    #   # generate phonon rattled structures
    #   fc2 = fcp.get_force_constants(supercell).get_fc_array(order=2, format='ase')
    #   structures_phonon_rattle = generate_phonon_rattled_structures(
    #       supercell, fc2, n_structures, T)
    #   write('structures_phonon_rattle_T{}.extxyz'.format(T), structures_phonon_rattle)  # noqa


# Module:
#   calculators

#   M: calculators
def get_aproximate_calculator(method='tblite'):
    """Returns requested calculator"""
    if method == 'ani':
        try:
            import torchani
        except Exception as e:
            print(type(e))
            print('ERROR: torchani is not installed. Bye for now!')
            sys.exit(1)
        calculator = torchani.models.ANI1ccx().ase()
    elif method == 'xtb':
        from xtb.ase.calculator import XTB
        calculator = XTB(method="GFN2-xTB")
    elif method == 'tblite':
        from tblite.ase import TBLite
        # calculator = TBLite(method="GFN1-xTB")
        calculator = TBLite(method="GFN2-xTB")

    return calculator


def preotimize(system, method='tblite', optimizer='bfgs', fmax=0.5, max_steps=None):
    """Performs a preoptimization of the system"""
    logging.info(f'Preoptimising with {method}; F_max = {fmax}; max steps = {max_steps}')
    optimize(system, method, optimizer, fmax, max_steps)


def optimize(system, method='tblite', optimizer='bfgs', fmax=0.01, max_steps=1000):
    """Performs an optimization of the system"""
    from ase.optimize import BFGS

    # system.pbc = np.array([False, False, False])
    system.calc = get_aproximate_calculator(method)

    # Geometry minimization of the system
    logging.info(f'Optimising with {method}; F_max = {fmax}; max steps = {max_steps}')
    if optimizer == 'bfgs':
        opt = BFGS(system, maxstep=max_steps)
    opt.run(fmax=fmax)


#   M: calculators
def sample_geometries(system, method='md'):
    if config:
        method = config.sampling_method

    if method == 'md':
        sampling_from_MD(system)
    if method == 'rattle':
        gen_rattled_geometries(system)

def sampling_from_MD(system, method='tblite', sampling_interval=20, temperature=500, md_steps=1000):
    """
    Performs a Molecular Dynamic using as calculator:
     - ani: torchani (ANAKIN-ME like Deep Learning potentials)
     - xtb: DFTB (does not work on periodic systems)
     - tblite: DFTB (works on periodic systems)
    It can preoptimize the geometry before start the MD
    """

    from ase.io import write
    from ase import units
    from ase.io.trajectory import Trajectory
    from ase.md.langevin import Langevin

    from math import floor

    if config:
        method = config.sampling_calculator
        sampling_interval = config.sampling_interval
        temperature = config.sampling_temperature
        md_steps = config.sampling_md_steps

    # system.pbc = np.array([False, False, False])
    system.calc = get_aproximate_calculator(method)

    def sample_geometry(format='extxyz'):
        """Samples a geometry from the MD"""
        import uuid
        calcfile = str(uuid.uuid4()).split('-')[4]

        # Define the subfolder path to save the sampled geometries
        sampling_subfolder = "MD_sampled_geometries"
        os.makedirs(sampling_subfolder, exist_ok=True)

        # Save the sampled geometry in the subfolder
        filepath = os.path.join(sampling_subfolder, calcfile + '.' + format)
        write(filepath, system)
        # write('trajectory.extxyz', system, append=True)

    dyn = Langevin(system, 1 * units.fs, temperature * units.kB, 0.2)

    traj = Trajectory('md.traj', 'w', system)
    dyn.attach(traj.write, interval=sampling_interval)
    dyn.attach(sample_geometry, interval=sampling_interval)

    logging.info(f'MD with {method} started:')
    logging.info(f'  T = {temperature}; MD steps = {md_steps}; sampling every {sampling_interval} steps')
    dyn.run(md_steps)
    sampled_structures = floor(md_steps/sampling_interval)
    logging.info(f'  MD finished. Saved in md.traj. Sampled {sampled_structures} strcutures')


#   M: calculators
def get_DFT_calculator_parameters():
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


#   M: calculators
def set_DFT_calculator_parameters():
    """ Set the calculator parameters from config object"""

    # TODO: When adding support for more calculators, come here and tweak it
    input_params, pseudos, kpts, command = get_DFT_calculator_parameters()
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


#   M: calculators
def get_QM_forces(path=None):
    """Calculates the DFT forces from the sampled geometries"""

    # TODO: Most of the lines of this function are to handle paths and filenames
    # Move to or use dirman module and leave here only the QM calculation part

    # If no path given, works on current folder
    cwd = os.getcwd()
    if path:
        os.chdir(path)

    # Has to be the same as the generated in the MD
    # Maybe create a parameter?
    sampling_subfolder = "MD_sampled_geometries"
    qm_forces_subfolder = "QM_forces"
    os.makedirs(qm_forces_subfolder, exist_ok=True)

    geometries = os.listdir(sampling_subfolder)

    # Go through the extxyz files in folder --> Read in ase format
    for geom in geometries:
        # work out the files and paths
        qm_force_filename = geom.replace('.extxyz', '.pwo')
        qm_force_filepath = os.path.join(qm_forces_subfolder, qm_force_filename)
        qm_input_filename = geom.replace('.extxyz', '.pwi')

        if os.path.exists(qm_force_filepath):
            continue  # geometry has been already calculated

        # Read the sampled geometry
        sampled_filepath = os.path.join(sampling_subfolder, geom)
        system = read(sampled_filepath)

        calculator = set_DFT_calculator_parameters()
        system.calc = calculator
        system.calc.write_input(system)

        os.rename('espresso.pwi', os.path.join(qm_forces_subfolder, qm_input_filename))

        # Calculate the forces of each geometry
        # TODO: Remove the if, and design a better workflow manipulation
        if (config.calculate_f):
            try:
                system.get_forces()
                logging.info(f'  Forces calculation finished: {qm_force_filepath}')
                # Save as .pwo and QM_filename.extxyz
                os.rename('espresso.pwo', qm_force_filepath)
                # FIXME: It does not write the .extxyz file of the DFT
                # extxyz_qm_file = qm_force_filepath.replace('.pwo', '.extxyz')
                # write(system, extxyz_qm_file)
            except:  # QE does not finish properly or optimization does not converge
                logging.error(f' QEspresso in {qm_force_filepath} finished with error.')

# Module:
#   visualize
#   M: visualize
def visualize(system):
    """Uses ASE GUI to visualize the system"""
    if config.visualize:
        view(system, repeat=config.vis_repeat)
