import numpy as np

# -------------------------
# Write files
# -------------------------
# write('system.in', system, format='espresso-in',
#       pseudopotentials={'C': '', 'O': ''})
# traj_file = 'geopt.traj'

# File manipulation stuff
#   M: utils
def get_parameter_from_pdb(pdb_file, pattern=None):
    with open(pdb_file, 'r') as file:
        for line in file:
            if pattern in line:
                return line

#   M: utils
def insert_line_in_pdb(pdb_file, new_line, insert_index):
    # Read the contents of the PDB file
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    # Insert the new line at the desired index
    lines.insert(insert_index, new_line + '\n')

    # Write the modified lines back to the PDB file
    with open(pdb_file, 'w') as file:
        file.writelines(lines)

#   M: utils
def get_xy_distance(a, b):
    """
    Calculates the distance between 2 vectors in their projections on
    the XY plane (ignores the Z coordinate)
    """
    return np.linalg.norm(a[0:2] - b[0:2])
