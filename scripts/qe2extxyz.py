import sys
import os
from ase.io import read, write

"""
Searchs for .pwo or else for .extxyz files and gathers the
geometries with their forces in one file geometry-forces.extxyz
"""

if len(sys.argv) < 2:
    command = 'find ./ -name "*pwo"  > tmp_out_files.txt'
    os.system(command)
    output_files = 'tmp_out_files.txt'
    if os.stat(output_files).st_size == 0:
        command = 'find ./ -name "*extxyz"  > tmp_out_files.txt'
        os.system(command)
        output_files = 'tmp_out_files.txt'
    if os.stat(output_files).st_size == 0:
        sys.exit("No relevant files found")
else:
    output_files = sys.argv[1]

with open(output_files) as f:
    for file in f:
        print('Reading file:', file[:-1])
        try:
            infile = read(file[:-1])
            outfile = write('geometry-forces.extxyz', infile, append=True)
        except Exception as e:
            print(f'Exception {type(e).__name__}: Error reading/writing file: {file}')

os.remove(output_files)
