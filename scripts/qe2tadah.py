import sys
import os
import numpy as np

from ase.io import read

"""
Searchs for .pwo or else for .extxyz files and gathers the
geometries with their forces and stresses in one tadah.dataset file
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

with open('tadah.dataset', 'a') as dataset:
    idx = 0
    with open(output_files) as f:
        print('Reading files...')
        for file in f:
            try:
                infile = read(file[:-1])
                idx += 1
            except Exception as e:
                print(f'Exception {type(e).__name__}: Error reading/writing file: {file}')

            #  Comment line: get File name to use as label
            label = file.split('/')
            dataset.write(label[-1][:-1] + f' {idx}\n')

            #  eweight fweight sweight: Should be given as parameters?
            eweight = 1.0
            fweight = 1.0
            sweight = 1.0
            dataset.write(f'{eweight} {fweight} {sweight}\n')

            #  Energy
            energy = infile.get_total_energy()
            dataset.write(f'{energy}\n')

            #  cell vector
            cell = infile.get_cell().array
            dataset.write(f'{cell[0][0]}   {cell[0][1]}    {cell[0][2]}\n')
            dataset.write(f'{cell[1][0]}   {cell[1][1]}    {cell[1][2]}\n')
            dataset.write(f'{cell[2][0]}   {cell[2][1]}    {cell[2][2]}\n')

            #  stress tensor
            try:
                stress = infile.get_stress(voigt=False)
                dataset.write(f'{stress[0][0]}   {stress[0][1]}    {stress[0][2]}\n')
                dataset.write(f'{stress[1][0]}   {stress[1][1]}    {stress[1][2]}\n')
                dataset.write(f'{stress[2][0]}   {stress[2][1]}    {stress[2][2]}\n')
            except Exception:
                dataset.write('0.0   0.0   0.0')
                dataset.write('0.0   0.0   0.0')
                dataset.write('0.0   0.0   0.0')

            #  Element px py pz fx fy fz
            chemical_symbols = infile.get_chemical_symbols()
            positions = infile.get_positions()
            try:
                forces = infile.get_forces()
            except Exception:
                forces = np.zeros((len(chemical_symbols), 3))

            for i in range(len(chemical_symbols)):
                ipositions = ' '.join(map(str, positions[i]))
                iforces = ' '.join(map(str, forces[i]))
                dataset.write(chemical_symbols[i] + " " + ipositions + "   " + iforces + '\n')

            # Blank line for separation between structures
            dataset.write('\n')
    print(f"Done! Read {idx} files")

if len(sys.argv) < 2:
    os.remove('tmp_out_files.txt')
