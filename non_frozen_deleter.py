'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023
'''
import MDAnalysis as mda
import numpy as np
import pandas as pd
import math
from typing import List, Tuple 
import Plane_fitting_functions as pff
import argparse

parser = argparse.ArgumentParser(description='find all molecules with distances between C3 and C14 < than a criterion')
parser.add_argument('-f', '--f', type = str, metavar = '', required = True, help = 'Trajectory file')
parser.add_argument('-s', '--s', type = str, metavar = '', required = True, help = 'Topology file')
parser.add_argument('-b', '--begin', type = int, metavar = '', required = False, help = 'Begin time in frame number', 
        default = 0)
parser.add_argument('-e', '--end', type = int, metavar = '', required = False, help = 'End time in frame number', 
        default = -2)
parser.add_argument('-st', '--stride', type = int, metavar = '', required = False, help = 'Stride in frame number', 
        default = 1)
parser.add_argument('-a', '--average', type = float, metavar = '', required = False, help = 'Predefined length between the two atoms in A. Default value is 14',
                    default = 14) 
parser.add_argument('-m', '--mobility', type = float, metavar = '', required = False, help = 'Criterion for mobility in A. Default value is 2',
                    default = 2) 
args = parser.parse_args()

def main(input_top: str, input_trr: str, begin: int, end: int, stride: int, average: float, mobility_criterion: float):
    no_file_extensnion = str(input_trr[:len(input_trr) - 4])

    u = mda.Universe(input_top, input_trr)

    sel_C8 = u.select_atoms("name C8").resids

    distances = []

    try:
        if end == -2 or end == len(u.trajectory):
            traj = u.trajectory[begin::stride]
        else:
            traj = u.trajectory[begin : end + 1: stride]
    except IndexError:
        print("index error\nexiting...")
        exit(7)

    for ts in traj:#u.trajectory:
        print(ts.time, end='\r')
        distances.append(pff.Head_Tail(u, 'C3', 'C14', avg_std = False))
    print('\n')
    distances = np.array(distances)

    # avg = np.average(distances)
    # std = np.std(distances)

    # print(avg, std, distances.shape)
    # print(len(distances))
    # print(len(distances[0]))

    not_frozen = []

    for i in range(len(distances[0])):
        for j in range(len(distances)):
            if distances[j][i] < average - mobility_criterion :
                not_frozen.append(sel_C8[i])
                break

    frozen = u.select_atoms("not resid {}".format(" ".join(str(_) for _ in not_frozen)))
    
    frozen_resids = frozen.select_atoms("name C8").resids 
    frozen.write(f'{no_file_extensnion}_frozen_only.gro')

    u.atoms = frozen

    with mda.Writer(f'{no_file_extensnion}_frozen_only.trr', frozen.n_atoms) as W:
        for ts in u.trajectory:
            print(ts.time, end='\r')
            W.write(frozen)

    with open('frozen-in-system.txt', 'a') as f:
        f.write(f'\n {no_file_extensnion} \n')
        f.write('frosen resids: \n')
        f.write(f'{frozen_resids} \n')
        f.write('number of frozen molecules' + '\n')
        f.write(f'{len(frozen_resids)} \n')
        f.write('disordered resids: \n')
        f.write(f'{not_frozen} \n')
        f.write('number of disordered molecules' + '\n')
        f.write(f'{len(not_frozen)} \n')
    
    # print(frozen.select_atoms("name C8").resids, len(frozen.select_atoms("name C8").resids))
    # print(not_frozen, len(not_frozen))

if __name__ == "__main__":
    main(input_top= args.s, input_trr= args.f, begin = args.begin, end = args.end, stride = args.stride, average=args.average, mobility_criterion=args.mobility)
