import MDAnalysis as mda 
import numpy as np
import pandas as pd
import math
from typing import List, Tuple 
import Plane_fitting_functions as pff



def main(input_top: str, input_trr: str):
    # input_top = str("Pure_HEX_cryst1_layer2\R2_pure_prod_300K-278K_1K_400ns-500ns_whole_cryst1_allone_layer2_Dict_allone_5d2.gro")
    # input_trr = str("Pure_HEX_cryst1_layer2\R2_pure_prod_300K-278K_1K_400ns-500ns_whole_cryst1_allone_layer2_Dict_allone_5d2.trr")
    no_file_extensnion = str(input_trr[:len(input_trr) - 4])
    u = mda.Universe(input_top, input_trr)

    num_col = int(len(u.select_atoms('resname HEX').residues.resids)*3)
    num_rows = len(u.trajectory)
    intercept_traj = np.empty((0, num_col))

    for ts in u.trajectory:
        sel_C8 = u.select_atoms("name C8").positions

        X = sel_C8.T[0].flatten()
        Y = sel_C8.T[1].flatten()
        Z = sel_C8.T[2].flatten()

        sel_C3 =  u.select_atoms("name C3").positions
        XC3 = sel_C3.T[0].flatten()
        YC3 = sel_C3.T[1].flatten()
        ZC3 = sel_C3.T[2].flatten()

        sel_C13 =  u.select_atoms("name C13").positions
        XC13 = sel_C13.T[0].flatten()
        YC13 = sel_C13.T[1].flatten()
        ZC13 = sel_C13.T[2].flatten()

        plane = pff.PCAPlane(sel_C8).plane()

        intercept_time = []
        for i in range(len(XC13)):
            point1 = np.array([XC3[i], YC3[i], ZC3[i]]) 
            point2 = np.array([XC13[i], YC13[i], ZC13[i]])

            # intercept = PointOfIntercept(point1, point2, plane).intercept()
            intercept = pff.PointOfIntercept(point1, point2, plane).rot_matrix()
            intercept = np.round(intercept, 2)
            intercept_time.append(intercept)
        intercept_time = np.array(intercept_time).reshape(-1, num_col)
        intercept_traj = np.append(intercept_traj, intercept_time, axis=0)    
    df = pd.DataFrame(intercept_traj)
    df.to_csv(f'{no_file_extensnion}_for_Voronoi.csv', index=True, header=False)
        
if __name__ == "__main__":
    main(input_top = '', input_trr = '')