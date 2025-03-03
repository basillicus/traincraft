import os
import logging
import numpy as np
from pathlib import Path

from ase.io import write
from ase.io import read

import calculeitors
from config import config


def sampling_from_MD(system, method='tblite', sampling_interval=20, temperature=500, timestep=1, md_steps=1000):
    """
    Performs a Molecular Dynamic using as calculator:
     - ani: torchani (ANAKIN-ME like Deep Learning potentials)
     - xtb: DFTB (does not work on periodic systems)
     - tblite: DFTB (works on periodic systems)
     - mace-off23 (elements: H,C,N,O,F,P,S,Cl,Br,and I; molecular or periodic systems)
     - mace-mp0 (elements 1-89, mainly for periodic systems)
    It can preoptimize the geometry before start the MD
    """

    from ase import units
    from ase.io.trajectory import Trajectory
    from ase.md.langevin import Langevin

    from math import floor

    if config:
        method = config.sampling_calculator
        sampling_interval = config.sampling_md_interval
        temperature = config.sampling_md_temperature
        timestep = config.sampling_md_timestep
        md_steps = config.sampling_md_steps

        # # FIXME: It has to be a nicer way to handle mace_params than this
        #
        # mace_params = {'mace_model_path': None,
        #                'device': None}
        # if method.startswith('mace'):
        #     mace_model_path = config.sampling_mace_model_path
        #     device = config.sampling_device
        #     mace_params = {'mace_model_path': mace_model_path,
        #                    'device': device}

        mace_params = {'mace_model_path': config.sampling_mace_model_path,
                       'device': config.sampling_device} if method.startswith('mace') else {}

    # system.pbc = np.array([False, False, False])
    system.calc = calculeitors.get_aproximate_calculator(method, **(mace_params or {}))

    def sample_geometry(format='extxyz'):
        """Samples a geometry from the MD"""
        # M: to dirman
        import uuid
        calcfile = str(uuid.uuid4()).split('-')[4]

        # Define the subfolder path to save the sampled geometries
        sampling_subfolder = "MD_sampled_geometries"
        os.makedirs(sampling_subfolder, exist_ok=True)

        # Save the sampled geometry in the subfolder
        filepath = os.path.join(sampling_subfolder, calcfile + '.' + format)
        write(filepath, system)
        # write('trajectory.extxyz', system, append=True)

    # M: to calculeitors?
    dyn = Langevin(system, timestep * units.fs, temperature * units.kB, 0.2)

    traj = Trajectory('md.traj', 'w', system)
    dyn.attach(traj.write, interval=sampling_interval)
    dyn.attach(sample_geometry, interval=sampling_interval)

    logging.info(f'MD with {method} started:')
    logging.info(f'  T = {temperature}; MD steps = {md_steps}; timestep  = {timestep} fs; sampling every {sampling_interval} steps')
    dyn.run(md_steps)
    sampled_structures = floor(md_steps/sampling_interval)
    logging.info(f'  MD finished. Saved in md.traj. Sampled {sampled_structures} strcutures')


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
    rattled_structures = config.sampling_rattle_nstructures
    rattle_std = config.sampling_rattle_std
    rattle_method = config.sampling_rattle_method

    # FIXME: It has to be a nicer way to handle mace_params than this
    mace_params = {'mace_model_path': None,
                   'device': None}
    if method.startswith('mace'):
        mace_model_path = config.sampling_mace_model_path
        device = config.sampling_device
        mace_params = {'mace_model_path': mace_model_path,
                       'device': device}

    calc = calculeitors.get_aproximate_calculator(method, **mace_params)

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
                write(f'structures_mc_rattle_{i}.extxyz', structures_mc_rattle)
            except Exception as e:
                logging.warn(e)
                pass

    rattle_goemetries = [f for f in os.listdir() if f.endswith('.extxyz')]

    if config.sampling_optimize_rattled:
        fmax = config.sampling_optimize_rattled_fmax
        max_steps = config.sampling_optimize_rattled_maxStep
        logging.info(f'Preoptimizing rattled geoemtries with Fmax: {fmax}; Max Steps: {max_steps}')
        for geom in rattle_goemetries:
            system = read(geom)
            logging.info(f"  Optimizing rattled structure:{geom} ")
            calculeitors._optimize(system, fmax=fmax, max_steps=max_steps)
            write(geom, system)

    # TODO: Try to get it from a centralize module: dirman
    sampling_subfolder = "rattle_sampled_geometries"
    p = Path(sampling_subfolder)
    os.makedirs(p, exist_ok=True)
    for f in rattle_goemetries:
        os.rename(f, p/f)
    #
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
