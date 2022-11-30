import MDAnalysis as mda 
import numpy as np
import pandas as pd
from typing import List, Tuple 
import Plane_fitting_functions as pff
from time import time

input_top = str("/home/stoyan/crystallites/R2-2Serfices-1K-400-500ns/0-500ns/R2-2Serfices-1K-400-500ns_whole_cryst1.gro")
input_trr = str("/home/stoyan/crystallites/R2-2Serfices-1K-400-500ns/0-500ns/R2-2Serfices-1K-0-500ns_cryst1.trr")
no_file_extensnion = str(input_trr[:len(input_trr) - 4])

start_t = 900
end_t = 3000

if input_trr == "None" or input_trr == "":
    u = mda.Universe(input_top)
else:
    u = mda.Universe(input_top, input_trr)

    avg_angl_nosmooth = []
    std_angl_nosmooth = []
    all_angl_nosmooth = []

    # stime = time()
   #pff.Progression_timer(0, len(u.trajectory[start_t:end_t])).progress_bar()
    for ts in u.trajectory[start_t:end_t]:

#       sel_C8 = u.select_atoms("name C8").positions

        plane = pff.PCAPlane(sel_C8).plane()

#       plane_normal = [0.6170, 0.3853, 0.6862]
        
        angl_nosmooth = pff.Vector_Plane_Angle(u, plane_normal)
        avg_angl_nosmooth.append(angl_nosmooth[0])
        std_angl_nosmooth.append(angl_nosmooth[1])
        all_angl_nosmooth.append(angl_nosmooth[2])
        ''' End Plane Fitting Without Smoothing '''
        
       #pff.Progression_timer(i+1, len(u.trajectory[start_t:end_t])).progress_bar()

    # Arrays with angles for the nonsmoothed plane fit
    avg_angl_nosmooth = np.array(avg_angl_nosmooth)
    std_angl_nosmooth = np.array(std_angl_nosmooth)

    all_angl_nosmooth = np.array(all_angl_nosmooth)
    all_angl_nosmooth_avg_time = np.average(all_angl_nosmooth, axis=0)
    all_angl_nosmooth_std_time = np.std(all_angl_nosmooth, axis=0)

    
    resids = mda.Universe(input_top).select_atoms("name C8").resids
    
    if end_t == '':
        end_t = len(u.trajectory)

    pd.DataFrame({'Tilt angle' : avg_angl_nosmooth, 
                'Tilt angle std' : std_angl_nosmooth}).to_csv(f'{no_file_extensnion}_avg_angle_over_time_PCA_{int(start_t/100)}-{int(end_t/100)}.csv', index=False)
    
    pd.DataFrame({'molecule number' : resids, 'Tilt angle' : all_angl_nosmooth_avg_time, 
                'Tilt angle std' : all_angl_nosmooth_std_time}).to_csv(f'{no_file_extensnion}_indiv_mol_angle_PCA_{int(start_t/100)}-{int(end_t/100)}.csv', index=False)
    # print(f'Elapsed Time: {((time() - stime)/60):.2f}')
