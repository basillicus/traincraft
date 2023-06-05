from tomlkit import load
from pathlib import Path


class Config:
    """
    A class for loading and accessing configuration parameters from a TOML file.

    Usage:
    config = Config('config.toml')
    config.read_parameters()
    value = config.get('key', default)

    Parameters:
    configfile (str): Path to the TOML configuration file.
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
        """Reads the required parameters.
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
        """ Called by read_parameters().
        Gets geometry parameters from config file"""

        self.geom_generation = self.config.get('geom_generation', 'packmol')

        self.n_structures = self.config['structures'].get('n_structures', 3)

        self.tolerance = self.config['molecules'].get('tolerance', 2.0)
        self.molec = self.config['molecules'].get('molecule', 'CO2')
        self.n_molecules = self.config['molecules'].get('n_molecules', 2)

        # Molecules manual manipulation
        self.compresion_factor = self.config['molecules'].get(
            'compresion_factor', 1.0)

        self.tilt_factor = self.config['molecules'].get('tilt_factor', None)
        self.molec_rot_x = self.config['molecules'].get('rot_x', None)
        self.molec_rot_y = self.config['molecules'].get('rot_y', None)
        self.molec_rot_z = self.config['molecules'].get('rot_z', None)

        self.molec_rot_axis = self.config['molecules'].get(
            'rot_axis', ['x', 'y', 'z'])

        # CNT
        self.cnt_n = self.config['cnt'].get('cnt_n', 8)
        self.cnt_m = self.config['cnt'].get('cnt_m', 0)
        self.cnt_l = self.config['cnt'].get('cnt_l', 2)
        self.cnt_bond = self.config['cnt'].get('cnt_bond', 1.43)
        self.cnt_gap = self.config['cnt'].get('cnt_gap', 4.0)
        self.cnt_constraints = self.config['cnt'].get('constraints', 'all')

    def _read_visualize_parameters(self):
        """ Called by read_parameters().
        Gets visualize parameters from config file"""
        self.visualize = self.config['visualize'].get('visualize', True)
        self.vis_repeat = self.config['visualize'].get('repeat', [2, 2, 2])

    def _read_calculator_parameters(self):
        """ Called by read_parameters().
        Gets calculator parameters from config file"""

        self.calculator = self.config['calculator'].get('calculator', 'qe')
        self.calculate_e = self.config['calculator'].get(
            'calculate_energies', False)
        self.calculate_f = self.config['calculator'].get(
            'calculate_forces', False)

        default_input_params_qe = {
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

        # define the pseudopotentials
        default_pseudos = {"C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                           "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
                           "N": "N.oncvpsp.upf",
                           "H": "H_ONCV_PBE-1.0.oncvpsp.upf",
                           }

        self.input_params = self.config['calculator'].get(
            'qe', default_input_params_qe)
        # print(self.input_params)
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
