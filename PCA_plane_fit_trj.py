import MDAnalysis as mda 
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
from typing import List, Tuple 
import Plane_fitting_functions as pff

def main(input_top: str, input_trr: str):
    input_top = str("R2-2Serfices-1K-400-500ns_whole_cryst6.gro")
    input_trr = str("R2-2Serfices-1K-400-500ns_whole_cryst6.gro")
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    u = mda.Universe(input_top, input_trr)

    avg_angl_nosmooth = []
    std_angl_nosmooth = []
    all_angl_nosmooth = []
    
    for ts in u.trajectory:

        sel_C8 = u.select_atoms("name C8").positions

        # x = sel_C8.T[0]
        # y = sel_C8.T[1]
        # z = sel_C8.T[2]

        plane = pff.PCAPlane(sel_C8).plane()

        Z_pred = []

        # for i in range(len(x)):
        #     prediction = (plane[0]*x[i] + plane[1]*y[i] - plane[3])/(-plane[2])
        #     Z_pred.append(prediction)

        # Z_pred = np.array(Z_pred)
        # plane_normal = [0.617, 0.385, 0.686]
        plane_normal = [plane[0], plane[1], plane[2]]
        
        angl_nosmooth = pff.Vector_Plane_Angle(u, plane_normal)
        avg_angl_nosmooth.append(angl_nosmooth[0])
        std_angl_nosmooth.append(angl_nosmooth[1])
        all_angl_nosmooth.append(angl_nosmooth[2])
        ''' End Plane Fitting Without Smoothing '''

    # Arrays with angles for the nonsmoothed plane fit
    avg_angl_nosmooth = np.array(avg_angl_nosmooth)
    std_angl_nosmooth = np.array(std_angl_nosmooth)

    all_angl_nosmooth = np.array(all_angl_nosmooth)
    all_angl_nosmooth_avg_time = np.average(all_angl_nosmooth, axis=0)
    all_angl_nosmooth_std_time = np.std(all_angl_nosmooth, axis=0)

    
    resids = mda.Universe(input_top).select_atoms("name C8").resids

    pd.DataFrame({'Tilt angle' : avg_angl_nosmooth, 
                'Tilt angle std' : std_angl_nosmooth}).to_csv(f'average_PCA_plane/{no_file_extensnion}_avg_angle_over_time_PCA_TEST.csv', index=False)
    
    pd.DataFrame({'molecule number' : resids, 'Tilt angle' : all_angl_nosmooth_avg_time, 
                'Tilt angle std' : all_angl_nosmooth_std_time}).to_csv(f'average_PCA_plane/{no_file_extensnion}_indiv_mol_angle_PCA_TEST.csv', index=False)

if __name__ == "__main__":
    main(input_top = '', input_trr = '')