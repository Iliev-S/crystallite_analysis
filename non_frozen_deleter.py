import MDAnalysis as mda
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import math
from typing import List, Tuple 
import Plane_fitting_functions as pff

def main(input_top: str, input_trr: str, mobility_criterion: float = 2):
    # input_top = str('R2-cryst4\R2-2Serfices-1K-400-500ns_whole_cryst4.gro')
    # input_trr = str('R2-cryst4\R2-2Serfices-1K-400-500ns_whole_cryst4_Dict_allone_7_crit.trr')
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    u = mda.Universe(input_top, input_trr)

    sel_C8 = u.select_atoms("name C8").resids

    distances = []

    for ts in u.trajectory:
        print(ts.time, end='\r')
        distances.append(pff.Head_Tail(u, 'C3', 'C14', avg_std = False))
    print('\n')
    distances = np.array(distances)

    avg = np.average(distances)
    std = np.std(distances)

    # print(avg, std, distances.shape)
    # print(len(distances))
    # print(len(distances[0]))

    not_frozen = []

    for i in range(len(distances[0])):
        for j in range(len(distances)):
            if distances[j][i] < avg - mobility_criterion :
                not_frozen.append(sel_C8[i])
                break

    frozen = u.select_atoms("not resid {}".format(" ".join(str(_) for _ in not_frozen)))

    frozen.write(f'{no_file_extensnion}_frozen_only.gro')

    u.atoms = frozen

    with mda.Writer(f'{no_file_extensnion}_frozen_only.trr', frozen.n_atoms) as W:
        for ts in u.trajectory:
            print(ts.time, end='\r')
            W.write(frozen)


    print(frozen.select_atoms("name C8").resids, len(frozen.select_atoms("name C8").resids))
    print(not_frozen, len(not_frozen))

if __name__ == "__main__":
    main(input_top = 'R2-289K/R2-278-289K_prod_whole_no_water_cryst1.gro', input_trr = 'R2-289K/R2-278-289K_prod_whole_no_water_cryst1.trr', mobility_criterion = 2)