'''If you are using this script, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023
'''
import MDAnalysis as mda
import numpy as np
from typing import List, Dict
import Plane_fitting_functions as pff
from MDAnalysis.transformations import rotateby
import argparse

parser = argparse.ArgumentParser(description='Imaging of molecules to reconstruct their real-space coordinates.')
parser.add_argument('-f', '--f', type = str, metavar = '', required = True, help = 'Trajectory file')
parser.add_argument('-s', '--s', type = str, metavar = '', required = True, help = 'Topology file')
parser.add_argument('-stm', '--straight_mol', type = str, metavar = '', required = True, help = 'Residue ID of a straight molecule (all trans configurations)')
parser.add_argument('-vac', '--vacuum', type = float, metavar = '', required = True, help = 'Vacuum width between the layers')
parser.add_argument('-ang', '--ang_corr', type = float, metavar = '', required = False, help = 'Add angle correction in degrees if needed. Default value is 0',
                    default = 0) 
parser.add_argument('-b', '--begin', type = int, metavar = '', required = False, help = 'Begin time in frame number', 
        default = 0)
parser.add_argument('-e', '--end', type = int, metavar = '', required = False, help = 'End time in frame number', 
        default = -2)
parser.add_argument('-st', '--stride', type = int, metavar = '', required = False, help = 'Stride in frame number', 
        default = 1)
args = parser.parse_args()

def main(input_top: str, input_trr: str, straight_mol: str, vacuum_width: float, begin: int, end: int, stride: int, ang_corr: float):
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    u = mda.Universe(input_top, input_trr)
    u2 = mda.Universe(input_top)
    
    id_to_index = pff.res_ind_dict(u)
    hex_len = 50
    surf_len = 65

    '''Rotation angle'''
    selC3 = u2.select_atoms(f"name C3 and resid {straight_mol}").positions[0]
    selC14 = u2.select_atoms(f"name C13 and resid {straight_mol}").positions[0]

    for_cross_prod = [selC14[0]-selC3[0], selC14[1]-selC3[1], selC14[2]-selC3[2]]

    Zvec = [0, 0, 1]
    cross_prod = np.cross(for_cross_prod, Zvec)
    cross_prod_norm = cross_prod/(np.sqrt(cross_prod[0]**2 + cross_prod[1]**2 + cross_prod[2]**2))
    rotation_angle = pff.Vector_Plane_Angle(u2, Zvec)[0]
    '''Rotation angle'''
    try:
        if end == -2 or end == len(u.trajectory):
            traj = u.trajectory[begin::stride]
        else:
            traj = u.trajectory[begin : end + 1: stride]
    except IndexError:
        print("index error\nexiting...")
        exit(7)
    with mda.Writer(f"{no_file_extensnion}_allone.trr", n_atoms = u.atoms.n_atoms, multiframe = True) as W:
        for ts in traj:#u.trajectory[begin : end + 1: stride]:
            print(ts.time, end='\r')
            box_X, box_Y, box_Z = u.dimensions[0], u.dimensions[1], u.dimensions[2]
            sel_all_pos = u.atoms.positions

            ag = u.atoms           
            rotated = rotateby(rotation_angle + ang_corr, cross_prod_norm, ag=ag)(ag)
            
            # rotated.atoms.write(f'probni/{no_file_extensnion}_frame{frame}.gro')
            resid_C8 = u.select_atoms("name C8").resids
            sel_C8 = []
            for i in resid_C8:
                center_of_mass = rotated.select_atoms(f'resid {i}').center_of_mass()
                sel_C8.append(center_of_mass[2])
            # Z_squish = sel_C8.T[2].flatten()


            rezen4eta = pff.Vacuum_Layer(sel_C8, resid_C8, vacuum_width).vacuum()
            

            n_layer = len(rezen4eta)
            
            main_body = pff.most_values(rezen4eta, n_layer)[0]
            main_body_index = list(rezen4eta.keys()).index(main_body)
            # print(rezen4eta.keys(), list(rezen4eta.keys()).index(main_body), main_body)
            # Positions check
            for ind, key in enumerate(rezen4eta.keys()):
                values = rezen4eta[key]
                pos_of_layer = main_body_index - list(rezen4eta.keys()).index(key) 
                #print(pos_of_layer)
                for i in values:
                    mol_name = u.select_atoms(f'resid {i} and name C8').resnames
                    if mol_name.size > 0 and mol_name == str('C16'):
                        mol_len = surf_len
                    elif mol_name.size > 0 and mol_name == str('HEX'): 
                        mol_len = hex_len
                    else: continue

                    counter = id_to_index[i]
                    for j in range(counter, counter + mol_len):
                        sel_all_pos[j][0] = sel_all_pos[j][0] + box_X*pos_of_layer
            u.atoms.positions = sel_all_pos
            # u.atoms.write(f'probni/{no_file_extensnion}_final_frame{frame}.gro')
            W.write(u.atoms)
    print('\n')

if __name__ == "__main__":
    main(input_top= args.s, input_trr= args.f, straight_mol= args.straight_mol, vacuum_width= args.vacuum, 
        ang_corr = args.ang_corr, begin = args.begin, end = args.end, stride = args.stride)
    print(f'File created: {str(args.s[:len(args.s) - 4])}_allone.trr')