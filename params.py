from tomlkit import load


class Config:
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

        being config and object of the Class Config:
            config = Config(configfile)
        """
        # ---------------------------
        # System Geometry parameters
        # ---------------------------
        self.molec = self.config['molecules'].get('molecule', 'CO2')
        self.n_molecules = self.config['molecules'].get('n_molecules', 2)
        self.add_tilt = self.config['molecules'].get('add_tilt', True)
        self.tilt_factor = self.config['molecules'].get('tilt_factor', 15)
        self.compresion_factor = self.config['molecules'].get(
            'compresion_factor', 1.0)

        # CNT
        self.cnt_n = self.config['cnt'].get('cnt_n', 8)
        self.cnt_m = self.config['cnt'].get('cnt_m', 0)
        self.cnt_l = self.config['cnt'].get('cnt_l', 2)
        self.cnt_gap = self.config['cnt'].get('cnt_gap', 4.0)

        # config
        self.random_seed = self.config['config'].get('random_seed', True)
        self.seed = self.config['config'].get('seed', 42)

        # ---------------------------
        # Calculation parameters
        # ---------------------------

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
                           "O": "O.pbe-n-kjpaw_psl.0.1.UPF"}

        self.input_params = self.config['calculator'].get(
            'qe', default_input_params_qe)
        # print(self.input_params)
        self.input_params.pop('pseudos')
        self.kpts = self.config['calculator']['qe'].get('kpts')
        # self.input_params.pop('kpts')
        self.pseudos = self.config['calculator']['qe'].get(
            'pseudos', default_pseudos)
