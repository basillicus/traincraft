import sys
import os
import logging
from params import Config

import random
import numpy as np

"""Initializes the program

Generates the config object and set up system stuff"""

default_config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inputs/default_config.toml')
configfile = default_config_file
if len(sys.argv) > 1:
    configfile = sys.argv[1]

# Create config object and get parameters from configfile
config = Config(configfile)
config.read_parameters()

# Initialize logging file
logging.basicConfig(filename=config.logfile, level=logging.INFO,
                    format='%(asctime)s [%(levelname)s] %(message)s')
# logging.captureWarnings(True)

def set_seed():
    """Duh!"""
    if config.random_seed:
        seed = random.randrange(2**32 - 2)
    else:
        seed = config.seed
    random.seed(seed)
    np.random.seed(seed)
    logging.info(f'seed: {seed}')


# Order the chaos
set_seed()
