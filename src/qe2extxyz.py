from ase.io import read, write

output_files = 'output_files.txt'

with open(output_files) as f:
    for file in f:
        print('Reading file:', file)
        try:
            infile = read(file[:-1])
            outfile = write('geometry-forces.extxyz', infile, append=True)
            print('Read file structure:', infile)
        except:
            print(f'Error reading/writing file: {file}')
