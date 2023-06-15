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

        print('  Reading geometry parameters...', end='')
        self._read_geometry_parameters()
        print('OK!')

        print('  Reading visualize parameters...', end='')
        self._read_visualize_parameters()
        print('OK!')

        print('  Reading calculator parameters...', end='')
        self._read_calculator_parameters()
        print('OK!')

        print('  Reading configuration parameters...', end='')
        self._read_configuration_parameters()
        print('OK!')

    def _read_geometry_parameters(self):
        """
        Called by read_parameters()

        Parameters
        ----------
        [structure]

        geom_generation : '*auto*' | 'manual'
            Determines how geometries will be generated
                - manual: controls the displacement and rotation of the
                  molecules with the keywords `rot_axis`, `rot_x`, `rot_y`, `rot_z`, etc
                - auto: will use Packmol to generate different geometries. Still
                  `rot_axis`, `rot_x`, `rot_y` and `rot_z`, can be used to provide Packmol
                  molecules with the desired orientation

        n_structures : *1*
            Number of different structures to be generated

        pbc : *True* | False | [False, False, True]
            Set the Periodic Boundary conditions *not yet implemented*

        Elements for geometry generation
        --------------------------------

        [molecules]
        molec = self.config['molecules'].get('molecule', 'CO2')
        n_molecules = self.config['molecules'].get('n_molecules', 2)
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
        cnt_l = [2]
        cnt_bond = [1.43]
        cnt_gap = [4.0]
        cnt_constraints = self.config['cnt'].get('constraints', 'all')
        """

        # Options: 'manual' | 'auto'
        self.geom_generation = self.config.get('geom_generation', 'auto')

        # Number of different structures to be generated
        self.n_structures = self.config['structures'].get('n_structures', 1)

        # TODO: Implement this keywords.
        #   if type is molecule_on_surface, or molecular, or molecular_mix, or laminar
        # Options: 'nanotube' | 'molecule' | ....
        self.structure_type = self.config['structures'].get('type', 'nanotube')
        # TODO: Implement this keyword.
        #  Allow for parcial periodic boundary conditions, or just molecular systems
        # Options: True | False |  [False, False, True]
        self.structure_pbc = self.config['structures'].get('pbc', True)

        # --------------------------------
        # ELEMENTS FOR GEOMETRY GENERATION
        # --------------------------------
        # Molecule
        # TODO:Allow molecules to be a list, so different molecules can be included
        # in the same system
        self.molec = self.config['molecules'].get('molecule', 'CO2')
        self.n_molecules = self.config['molecules'].get('n_molecules', 2)
        self.tolerance = self.config['molecules'].get('tolerance', 2.0)

        # Molecules manual manipulation
        self.compresion_factor = self.config['molecules'].get(
            'compresion_factor', 1.0)

        self.molec_rot_x = self.config['molecules'].get('rot_x', None)
        self.molec_rot_y = self.config['molecules'].get('rot_y', None)
        self.molec_rot_z = self.config['molecules'].get('rot_z', None)

        self.molec_rot_axis = self.config['molecules'].get('rot_axis', ['x', 'y', 'z'])

        # CNT
        self.cnt_n = self.config['cnt'].get('cnt_n', 8)
        self.cnt_m = self.config['cnt'].get('cnt_m', 0)
        self.cnt_l = self.config['cnt'].get('cnt_l', 2)
        self.cnt_bond = self.config['cnt'].get('cnt_bond', 1.43)
        self.cnt_gap = self.config['cnt'].get('cnt_gap', 4.0)
        self.cnt_constraints = self.config['cnt'].get('constraints', 'all')

        # How sample geometries: 'md' | 'rattle'
        self.sampling_method = self.config['sampling'].get('method', 'md')
        self.sampling_calculator = self.config['sampling'].get('calculator', 'tblite')

        # sampling.md:
        self.sampling_interval = self.config['sampling']['md'].get('sampling_interval', 20)
        self.sampling_temperature = self.config['sampling']['md'].get('temperature', 500)
        self.sampling_md_steps = self.config['sampling']['md'].get('md_steps', 1000)

        # sampling.rattle
        self.rattle_method = self.config['sampling']['rattle'].get('method', 'mc')  # 'mc' | 'standard'
        self.rattle_structures = self.config['sampling']['rattle'].get('rattle_structures', 10)
        self.rattle_std = self.config['sampling']['rattle'].get('rattle_std', 0.12)

    def _read_visualize_parameters(self):
        """ Called by read_parameters().
        Gets visualize parameters from config file"""
        self.visualize = self.config['visualize'].get('visualize', False)
        self.vis_repeat = self.config['visualize'].get('repeat', [2, 2, 2])

    def _read_calculator_parameters(self):
        """ Called by read_parameters().
        Gets calculator parameters from config file"""

        self.calculator = self.config['calculator'].get('calculator', 'qe')
        # TODO: Remove calcualte_e
        # self.calculate_e = self.config['calculator'].get('calculate_energies', False)
        self.calculate_f = self.config['calculator'].get('calculate_forces', False)

        default_input_params_qe = {
            "ecutwfc": 45,     # plane-wave wave-function cutoff
            "ecutrho": 180,    # density wave-function cutoff,
            "conv_thr": 1e-6,  # DFT self-consistency convergence
            "pseudo_dir": "~/pseudos/qe/SSSP_1.1.2_PBE_precision/",
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
        self.input_params.pop('pseudos')
        self.kpts = self.config['calculator']['qe'].get('kpts', None)
        # self.input_params.pop('kpts')
        self.pseudos = self.config['calculator']['qe'].get(
            'pseudos', default_pseudos)

    def _read_configuration_parameters(self):
        """ Called by read_parameters().
        Gets configuration parameters from config file"""

        # config
        self.random_seed = self.config['config'].get('random_seed', True)
        self.seed = self.config['config'].get('seed', 42)

        # Paths
        self.prefix = self.config['config'].get('folder_prefix', 'cnt')
        self.cwd = Path.cwd()

        outdir_str = self.config['config'].get('outdir')
        if outdir_str is not None:
            self.outdir = Path(outdir_str)
        else:
            self.outdir = self.cwd
        self.logfile = self.config['config'].get('logfile', 'logfile.log')

        # PERF: nproc is still experimental and still needs to be tested in
        # different systems with different environments

        # Parallelization
        self.nproc = self.config['config'].get('nproc', None)
