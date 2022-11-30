from ast import arg
import MDAnalysis as mda
import numpy as np 
import pandas as pd
import Plane_fitting_functions as pff
import sys

input_top = str('orthorhombic/HEX_R1_no_rot_prod_250-500ns_whole.gro')
no_file_extensnion = str(input_top[:len(input_top) - 4])

u = mda.Universe(input_top)

with open('orthorhombic/no_rot_selections.txt', 'r') as f:
    lines = f.readlines()
    layer = 1
    with mda.selections.gromacs.SelectionWriter(f'{no_file_extensnion}.ndx', mode='w') as ndx:
        for line in lines:
            sel_layer = u.select_atoms(f'{line}')
            ndx.write(sel_layer, name=f'layer{layer}')
            layer += 1