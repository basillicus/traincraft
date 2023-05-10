from ase.db import connect
from ase.io import read
# from glob import glob

path_to_database = '/home/david/calc/cnt/cnt_CO2.db'
file_with_output_files = 'output_files.txt'
path_to_output_files = []

db = connect(path_to_database)
with open(file_with_output_files) as f:
    for line in f:
        atoms = read(line[:-1])  # remove the \n character
        # db.write(atoms, filename=fname, format='Espresso')
        db.write(atoms, filename=line, format='Espresso')

# for fname in glob(path_to_output_files):
#     atoms = read(fname)
#     db.write(atoms, filename=fname, format='Espresso')
