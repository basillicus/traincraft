import sys
from params import Config

configfile = 'cnt_generator.toml'
if len(sys.argv) > 1:
    configfile = sys.argv[1]

config = Config(configfile)
config.read_parameters()
