"""All stuff related to calculators and calculations"""
import sys
import os
import logging

from ase.io import read
from ase.io import write
from ase.calculators.espresso import Espresso

from config import config
import gengeom

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


def preotimize(system, method='tblite', fmax=0.5, max_steps=None, optimizer='bfgs', trajectory='preopt.traj'):
    """Performs a preoptimization of the system"""
    if config:
        method = config.preopt_calculator
        fmax = config.preopt_fmax
        max_steps = config.preopt_maxSteps

    logging.info(f'Preoptimising with {method}; F_max = {fmax}; max steps = {max_steps}')
    _optimize(system, method=method, fmax=fmax, max_steps=max_steps, optimizer=optimizer)


def _optimize(system, method='tblite', fmax=0.01, max_steps=None, optimizer='bfgs', trajectory='optim.traj'):
    """Performs an optimization of the system"""
    from ase.optimize import BFGS

    # system.pbc = np.array([False, False, False])
    system.calc = get_aproximate_calculator(method)

    # Geometry minimization of the system
    if optimizer == 'bfgs':
        # opt = BFGS(system, maxstep=max_steps)
        opt = BFGS(system, maxstep=max_steps, trajectory=trajectory)
    opt.run(fmax=fmax)


#   M: calculators
def get_DFT_calculator_parameters():
    """ Get the calculator parameters from config object"""

    # TODO: When adding support for more calculators, come here and tweak it
    if config.calculator == 'qe':
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
    if config.calculator == 'otherCalculators':
        pass


#   M: calculators
def set_DFT_calculator_parameters():
    """ Set the calculator parameters from config object"""

    # TODO: When adding support for more calculators, come here and tweak it
    if config.calculator == 'qe':
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
    if config.calculator == 'otherCalculators':
        pass

    return calc


#   M: calculators
def get_QM_forces_from_sampled(path=None, sampled_by=None):
    """Calculates the DFT forces from the sampled geometries"""

    # If no path given, works on current folder
    cwd = os.getcwd()
    if path:
        os.chdir(path)
    if path is None:
        sampling_subfolder = "."

    # Has to be the same as the generated in the MD
    # Maybe create a parameter?
    if sampled_by == 'md':
        sampling_subfolder = "MD_sampled_geometries"
    elif sampled_by == 'rattle':
        sampling_subfolder = "rattle_sampled_geometries"

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

        # From this, once we have system we could call get_QM_forces(system)
        # but the logging and file handling may require some care
        calculator = set_DFT_calculator_parameters()
        system.calc = calculator
        system.calc.write_input(system)

        os.rename('espresso.pwi', os.path.join(qm_forces_subfolder, qm_input_filename))

        # Calculate the forces of each geometry
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

def get_QM_forces(system):
    """Calculates the forces of the given system

    Parameters
    ----------
    system: (ASE Atoms)
    """
    cwd = os.getcwd()
    calculator = set_DFT_calculator_parameters()
    system.calc = calculator
    # it may happen that it writes it twice
    system.calc.write_input(system)
    # Calculate the forces of each geometry
    try:
        logging.info(f'  Calculation of  forces started with {config.calculator}')
        system.get_forces()
        logging.info('  Forces calculation finished')
        # Save as .pwo and QM_filename.extxyz
        # os.rename('espresso.pwo', qm_force_filepath)
        # FIXME: It does not write the .extxyz file of the DFT
        # extxyz_qm_file = qm_force_filepath.replace('.pwo', '.extxyz')
        # write(system, extxyz_qm_file)
    except:  # QE does not finish properly or optimization does not converge
        logging.error(f' QEspresso in {cwd} finished with error.')
