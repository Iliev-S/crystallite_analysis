'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023
'''
import MDAnalysis as mda 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Tuple 
import Plane_fitting_functions as pff
import random
import argparse

parser = argparse.ArgumentParser(description='Fits a plane to the coordinates of the C8 atoms in a trajectory file and measures the angle of the molecules with respect to the plane. \n Outputs two files: One with the average tilt angle evolving in time and one with the tilt angle of each individual molecule averaged over time.')
parser.add_argument('-f', '--f', type = str, metavar = '', required = True, help = 'Trajectory file')
parser.add_argument('-s', '--s', type = str, metavar = '', required = True, help = 'Topology file')
parser.add_argument('-b', '--begin', type = int, metavar = '', required = False, help = 'Begin time in frame number', 
        default = 0)
parser.add_argument('-e', '--end', type = int, metavar = '', required = False, help = 'End time in frame number', 
        default = -2)
parser.add_argument('-st', '--stride', type = int, metavar = '', required = False, help = 'Stride in frame number', 
        default = 1)
parser.add_argument('-sp', '--saveplane', type = bool, metavar = '', required = False, help = 'Fits a plane to the coordinates from the provided .gro (coordinates) file and saves it to a .csv file.',
                    default = False) 
args = parser.parse_args()

def main(input_top: str, input_trr: str, begin: int, end: int, stride: int, save_plane: bool = False):
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    if save_plane is True:
        u2 = mda.Universe(input_top)
        
        sel_C8 = u2.select_atoms("name C8").positions

        x = sel_C8.T[0]
        y = sel_C8.T[1]
        z = sel_C8.T[2]
        
        min_x, max_x = min(x), max(x)
        min_y, max_y = min(y), max(y)

        plane = pff.PCAPlane(sel_C8).plane()

        Z_pred = []
        range_x = []
        range_y = []

        for i in range(1, 30001):
            ran_x = random.uniform(min_x, max_x)
            ran_y = random.uniform(min_y, max_y)
            prediction = (plane[0]*ran_x + plane[1]*ran_y - plane[3])/(-plane[2])
            Z_pred.append(prediction)
            range_x.append(ran_x)
            range_y.append(ran_y)

        Z_pred = np.array(Z_pred)
        range_x = np.array(range_x)
        range_y = np.array(range_y)

        pd.DataFrame({'X' : range_x, 'Y' : range_y, 'Z predicted' : Z_pred}).to_csv(f'{no_file_extensnion}_plane_visualisation.csv', index=False)
    
    
    u = mda.Universe(input_top, input_trr)

    avg_angl_nosmooth = []
    std_angl_nosmooth = []
    all_angl_nosmooth = []
    
    try:
        if end == -2 or end == len(u.trajectory):
            traj = u.trajectory[begin::stride]
        else:
            traj = u.trajectory[begin : end + 1: stride]
    except IndexError:
        print("index error\nexiting...")
        exit(7)

    for ts in traj:

        sel_C8 = u.select_atoms("name C8").positions
        
        plane = pff.PCAPlane(sel_C8).plane()

        Z_pred = []

        plane_normal = [plane[0], plane[1], plane[2]]
        
        angl_nosmooth = pff.Vector_Plane_Angle(u, plane_normal)
        avg_angl_nosmooth.append(angl_nosmooth[0])
        std_angl_nosmooth.append(angl_nosmooth[1])
        all_angl_nosmooth.append(angl_nosmooth[2])

    # Arrays with angles for the nonsmoothed plane fit
    avg_angl_nosmooth = np.array(avg_angl_nosmooth)
    std_angl_nosmooth = np.array(std_angl_nosmooth)

    all_angl_nosmooth = np.array(all_angl_nosmooth)
    all_angl_nosmooth_avg_time = np.average(all_angl_nosmooth, axis=0)
    all_angl_nosmooth_std_time = np.std(all_angl_nosmooth, axis=0)

    
    resids = mda.Universe(input_top).select_atoms("name C8").resids

    pd.DataFrame({'Tilt angle' : avg_angl_nosmooth, 
                'Tilt angle std' : std_angl_nosmooth}).to_csv(f'{no_file_extensnion}_avg_angle_over_time.csv', index=False)
    
    pd.DataFrame({'molecule number' : resids, 'Tilt angle' : all_angl_nosmooth_avg_time, 
                'Tilt angle std' : all_angl_nosmooth_std_time}).to_csv(f'{no_file_extensnion}_indiv_mol_angle.csv', index=False)

if __name__ == "__main__":
    main(input_top= args.s, input_trr= args.f, begin = args.begin, end = args.end, stride = args.stride, save_plane = args.saveplane)