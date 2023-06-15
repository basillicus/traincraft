import sys
import logging
from params import Config

"""Initializes the program

Generates the config object and set up loggin system"""

configfile = 'cnt_generator.toml'
if len(sys.argv) > 1:
    configfile = sys.argv[1]

# Create config object and get parameters from configfile
config = Config(configfile)
config.read_parameters()

# Initialize logging file
logging.basicConfig(filename=config.logfile, level=logging.INFO,
                    format='%(asctime)s [%(levelname)s] %(message)s')
