import sys

from tomlkit import load
from pathlib import Path


class Config:
    """
    Loading and accessing configuration parameters from a TOML file.

    Parameters:
    ----------

    configfile (str): Path to the TOML configuration file.

    Usage
    -----

    config = Config('config.toml')
    config.read_parameters()
    value = config.get('key', default)
    """

    def __init__(self, configfile):
        self.configfile = configfile
        self._config = None

    @property
    def config(self):
        if self._config is None:
            self._config = self._load_config()
        return self._config

    def _load_config(self):
        with open(self.configfile, 'r') as f:
            return load(f)

    def get(self, key, default=None):
        return self.config.get(key, default)

    def read_parameters(self):
        """Reads the required parameters from configfile

        Needs to be called after initializing a Config class.
        Other parameters should be able to be read when required with the
        method:
            config.get(key='key1.key2', default )

        where config is an object of the Class Config:
            config = Config(configfile)
        """

        self.title = self.config.get('title', 'untitled')
        print('Reading parameters for work:', self.title)

        # TODO: Implement parameter optimizaztion
        # self._read_DFT_parameters_optimization()

        print('  Reading geometry parameters... ', end='')
        self._read_geometry_parameters()
        print('OK!')

        print('  Reading visualize parameters... ', end='')
        self._read_visualize_parameters()
        print('OK!')

        print('  Reading calculator parameters... ', end='')
        self._read_calculator_parameters()
        print('OK!')

        print('  Reading configuration parameters... ', end='')
        self._read_configuration_parameters()
        print('OK!')

    def _read_geometry_parameters(self):
        """
        Called by read_parameters()

        Parameters
        ----------
        [structure]

        geom_generation : (str) or (bool) '*auto*' | 'manual' | False
            Determines how geometries will be generated
                - manual: controls the displacement and rotation of the
                  molecules with the keywords `rot_axis`, `rot_x`, `rot_y`, `rot_z`, etc
                - auto: will use Packmol to generate different geometries. Still
                  `rot_axis`, `rot_x`, `rot_y` and `rot_z`, can be used to provide Packmol
                  molecules with the desired orientation
                - False: Do not generate geometries

        n_structures : int
            Number of different structures to be generated

        pbc : *True* | False | [False, False, True]
            Set the Periodic Boundary conditions *not yet implemented*

        Elements for geometry generation
        --------------------------------

        [molecules]
        molec = self.config['molecules'].get('molecule', 'CO2')
        n_molecules = self.config['molecules'].get('n_molecules', 1)
        tolerance = self.config['molecules'].get('tolerance', 2.0)

        # Molecules manual manipulation
        compresion_factor = self.config['molecules'].get(
            'compresion_factor', 1.0)

        molec_rot_x = self.config['molecules'].get('rot_x', None)
        molec_rot_y = self.config['molecules'].get('rot_y', None)
        molec_rot_z = self.config['molecules'].get('rot_z', None)

        molec_rot_axis = self.config['molecules'].get(
            'rot_axis', ['x', 'y', 'z'])

        [cnt]
        cnt_n = [8]
        cnt_m = [0]
        cnt_l = [1]
        cnt_bond = [1.44]
        cnt_gap = [4.0]
        cnt_constraints = self.config['cnt'].get('constraints', 'all')

        [crystal]
        ...

        """

        # Options: 'manual' | 'auto'
        self.geom_generation = self.config.get('geom_generation', 'auto')

        # Number of different structures to be generated
        self.n_structures = self.config['structures'].get('n_structures', 1)

        # TODO: Implement this keywords.
        #   if type is molecule_on_surface, or molecular, or molecular_mix, or laminar
        # Options: 'nanotube' | 'molecule' | 'crystal' | 'surface' | 'mol on surf' ...
        self.structure_type = self.config['structures'].get('type', 'nanotube')
        # TODO: Implement this keyword.
        #  Allow for parcial periodic boundary conditions, or just molecular systems
        # Options: True | False |  [False, False, True]
        self.structure_pbc = self.config['structures'].get('pbc', True)

        # --------------------------------
        # ELEMENTS FOR GEOMETRY GENERATION
        # --------------------------------

        # MOLECULE
        self.molecules = self.config.get('molecules')
        if self.molecules is not None:
            # TODO:Allow molec to be a list, so different molecules
            # can be included in the same system
            # Allow to give a input file with a geometry
            self.molec = self.config['molecules'].get('molecule', None)
            if self.molec is None:
                print("ERROR: key 'molecule' in table [molecules] not specified")
                print("Some examples: H2O, CO, CO2, methane... you choose")
                sys.exit(1)

            self.n_molecules = self.config['molecules'].get('n_molecules', 1)
            self.tolerance = self.config['molecules'].get('tolerance', 2.0)

            # Molecules manual manipulation
            self.compresion_factor = self.config['molecules'].get(
                'compresion_factor', 1.0)

            self.molec_rot_x = self.config['molecules'].get('rot_x', None)
            self.molec_rot_y = self.config['molecules'].get('rot_y', None)
            self.molec_rot_z = self.config['molecules'].get('rot_z', None)

            self.molec_rot_axis = self.config['molecules'].get('rot_axis', ['x', 'y', 'z'])

        # CNT
        self.cnt = self.config.get('cnt')
        if self.cnt is not None:
            self.cnt_n = self.config['cnt'].get('cnt_n', 8)
            self.cnt_m = self.config['cnt'].get('cnt_m', 0)
            self.cnt_l = self.config['cnt'].get('cnt_l', 1)
            self.cnt_bond = self.config['cnt'].get('cnt_bond', 1.44)
            self.cnt_gap = self.config['cnt'].get('cnt_gap', 4.0)
            self.cnt_constraints = self.config['cnt'].get('constraints', None)

        # SAMPLING
        # How sample geometries: 'md' | 'rattle'
        # By default do not do sampling
        self.do_sampling = False
        if self.config.get('sampling'):
            # If the Table is present, check for the flag, if not present set it to True
            self.do_sampling = self.config['sampling'].get('do_sampling', True)
            # Options: 'md' | 'rattle'
            self.sampling_method = self.config['sampling'].get('method', 'md')
            # opts: tblite | xtb | ani
            self.sampling_calculator = self.config['sampling'].get('calculator', 'tblite')
            self.sampling_device = self.config['sampling'].get('device', 'cpu')
            self.sampling_mace_model_path = self.config['sampling'].get('mace_model_path', None)

            # sampling.md:
            if self.config['sampling'].get('md'):
                self.sampling_md_interval = self.config['sampling']['md'].get('sampling_interval', 20)
                self.sampling_md_temperature = self.config['sampling']['md'].get('temperature', 500)
                self.sampling_md_timestep = self.config['sampling']['md'].get('timestep', 0.5)
                self.sampling_md_steps = self.config['sampling']['md'].get('md_steps', 1000)

            # sampling.rattle
            # opts: 'mc' | 'standard' : Monte Carlo or Standar (chech HiPhive documentation for more info)
            if self.config['sampling'].get('rattle'):
                self.sampling_rattle_method = self.config['sampling']['rattle'].get('method', 'mc')
                self.sampling_rattle_nstructures = self.config['sampling']['rattle'].get('rattle_structures', 10)
                self.sampling_rattle_std = self.config['sampling']['rattle'].get('rattle_std', 0.12)

                self.sampling_optimize_rattled = False
                # TODO: Allow option to optimize rattled strucutres
                if self.config['sampling']['rattle'].get('optimize'):
                    self.sampling_optimize_rattled = True
                    self.sampling_optimize_rattled_fmax = self.config['sampling']['rattle']['optimize'].get('fmax', 1.0)
                    self.sampling_optimize_rattled_maxStep = self.config['sampling']['rattle']['optimize'].get('max_step', 0.1)

    def _read_visualize_parameters(self):
        # How sample geometries: 'md' | 'rattle'
        """ Called by read_parameters().
        Gets visualize parameters from config file"""
        self.visualize = self.config.get('visualize', False)
        if self.visualize:
            self.visualize = self.config['visualize'].get('visualize', False)
            self.vis_repeat = self.config['visualize'].get('repeat', [2, 2, 2])

    def _read_calculator_parameters(self):
        """ Called by read_parameters().
        Gets calculator parameters from config file"""

        # DFT Calculator parameters
        self.calculator = self.config.get('calculator')
        if self.calculator is not None:
            self.calculate_f = self.config['calculator'].get('calculate_forces', True)
            self.calculator = self.config['calculator'].get('calculator', 'qe')
            self.calculator_writeInput = self.config['calculator'].get('write_input_file', True)

        if self.calculator == 'qe':
            default_input_params_qe = {
                "ecutwfc": 50,     # plane-wave wave-function cutoff
                "ecutrho": 400,    # density wave-function cutoff,
                "conv_thr": 1e-6,  # DFT self-consistency convergence
                "pseudo_dir": "/home/david/pseudos/qe/SSSP_1.1.2_PBE_precision/",
                "vdw_corr": "xdm",
                "occupations": 'smearing',   # Add smearing
                "smearing": 'cold',  # smearing kind
                "degauss": 0.02,     # smearing amount
                "tprnfor": True,     # Print forces
                "tstress": True      # Print stress tensor
            }

            # define the pseudopotentials
            default_pseudos = {"C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                               "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
                               "N": "N.oncvpsp.upf",
                               "H": "H_ONCV_PBE-1.0.oncvpsp.upf",
                               }

            self.input_params = self.config['calculator'].get('qe', default_input_params_qe)
            # self.input_params.pop('pseudos')
            try:
                self.kpts = self.config['calculator']['qe'].get('kpts', None)
            except Exception as e:
                self.kpts = None
                pass
            # self.input_params.pop('kpts')
            try:
                self.pseudos = self.config['calculator']['qe'].get(
                    'pseudos', default_pseudos)
            except Exception as e:
                self.pseudos = default_pseudos

        # Preoptimize parameters
        self.do_preoptimize = False
        if self.config.get('preoptimize'):
            self.do_preoptimize = self.config['preoptimize'].get('do_preoptimize', True)
            # opt: tblite | xtb | mace-off23 | mace-mp0 | ani
            self.preopt_calculator = self.config['preoptimize'].get('calculator', 'tblite')
            self.preopt_fmax = self.config['preoptimize'].get('fmax', 2.0)
            self.preopt_maxSteps = self.config['preoptimize'].get('max_step', 0.25)
            self.preopt_mace_model_path = self.config['preoptimize'].get('mace_model_path', None)

    def _read_configuration_parameters(self):
        """ Called by read_parameters().
        Gets configuration parameters from config file"""

        # This variables need to exist, that is why we asign them a defaulf value
        # regardless of the corresponding table is present in the TOML file or not
        self.cwd = Path.cwd()
        self.outdir = self.cwd
        self.prefix = ''
        self.random_seed = True
        self.logfile = 'logfile.log'

        # Now check if the table [config] exists in the toml file.
        # If it exists, the default values may (or may not) be override if given
        # in the config file (or not)
        self.configuration = self.config.get('config')
        if self.configuration is not None:
            # config
            self.random_seed = self.config['config'].get('random_seed', True)
            self.seed = self.config['config'].get('seed', 42)

            # Paths
            outdir_str = self.config['config'].get('outdir', self.cwd)
            self.outdir = Path(outdir_str)
            self.prefix = self.config['config'].get('folder_prefix', '')
            self.logfile = self.config['config'].get('logfile', 'logfile.log')

            # PERF: nproc is still experimental. Not needed for QE >= 7
            # in versions of QE < 7 set nproc to the number of cores availables

            # Parallelization
            self.nproc = self.config['config'].get('nproc', None)
