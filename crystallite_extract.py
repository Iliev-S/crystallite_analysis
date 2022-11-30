from cgi import print_directory
import MDAnalysis as mda 
import numpy as np
from numpy.core.arrayprint import _extendLine
from numpy.core.fromnumeric import shape
import pandas as pd
from pandas.core.algorithms import mode
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import Plane_fitting_functions


def main(input_top: str, input_trr: str, cryst_sel: str):
    '''input_top: Name of the gro file
       input_trr: Name of the trr file
       cryst_sel: Name of a text file where each row cintains the resids of a crystallite.'''
    # input_top = 'R2-2Serfices-1K-400-500ns_whole.gro'
    # input_trr = 'R2-2Serfices-1K-400-500ns_whole.trr'
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    u = mda.Universe(input_top, input_trr)

    with open(cryst_sel, 'r') as f:
        lines = f.readlines()
        cryst_num = 1

        for line in lines:
            cryst_resids = line.split()
            print(cryst_resids)
            init_len = len(cryst_resids)
            for element in range(init_len):
                name = u.select_atoms(f'resid {cryst_resids[element]} and name C8').resnames
                if name.size > 0 and name == str('C16'):
                    eo2 = int(cryst_resids[element]) + 1
                    cryst_resids.append(eo2)
                
            cryst = u.select_atoms("resid {}".format(" ".join(str(_) for _ in cryst_resids)))
            cryst.write(f'{no_file_extensnion}_cryst{cryst_num}.gro')
            # cryst_num += 1
            with mda.Writer(f'{no_file_extensnion}_cryst{cryst_num}.trr', n_atoms = cryst.n_atoms, multiframe = True) as W:
                for ts in u.trajectory:
                    W.write(cryst)

                cryst_num += 1

if __name__ == "__main__":
    main(input_top= '', input_trr= '', cryst_sel= '')