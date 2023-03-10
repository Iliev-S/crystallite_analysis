'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023, 5, 2327. doi: 10.3390/molecules28052327
'''
import MDAnalysis as mda
import numpy as np 
from MDAnalysis.transformations import rotateby
import Plane_fitting_functions as pff
import argparse

parser = argparse.ArgumentParser(description='Fits a plane to the coordinates of the C8 atoms in a trajectory file and measures the angle of the molecules with respect to the plane. \n Outputs two files: One with the average tilt angle evolving in time and one with the tilt angle of each individual molecule averaged over time.')
parser.add_argument('-f', '--f', type = str, metavar = '', required = False, help = 'Trajectory file', default= '')
parser.add_argument('-s', '--s', type = str, metavar = '', required = True, help = 'Topology file')
parser.add_argument('-b', '--begin', type = int, metavar = '', required = False, help = 'Begin time in frame number', 
        default = 0)
parser.add_argument('-e', '--end', type = int, metavar = '', required = False, help = 'End time in frame number', 
        default = -2)
parser.add_argument('-st', '--stride', type = int, metavar = '', required = False, help = 'Stride in frame number', 
        default = 1)
parser.add_argument('-ar', '--axis_of_rotation', type = str, metavar = '', choices=['x', 'y', 'z'], required = False, help = 'Axis to rotate about. Choice of x, y, or z. Default is z',
                     default='z') 
parser.add_argument('-rd', '--rot_dir', type = int, metavar = '', choices=[1, -1], required = False, help = 'Clockwise/Anticlockwise rotation. 1 - Clockwise, -1 - Anticlockwise. Default is 1', 
        default = 1)
parser.add_argument('-stm', '--straight_mol', type = str, metavar = '', required = True, help = 'Residue ID of a straight molecule (all trans configurations)')
args = parser.parse_args()

def main(input_top: str, input_trr: str, begin: int, end: int, stride: int, straight_mol: str, axis_of_rotation: str, rot_dir: int):
    no_file_extensnion = str(input_top[:len(input_top) - 4])

    if axis_of_rotation == 'x':
        Zvec = [1, 0, 0]
    elif axis_of_rotation == 'y':
        Zvec = [0, 1, 0]
    else:
        Zvec = [0, 0, 1]

    if input_trr == 'None' or input_trr == '':
        u = mda.Universe(input_top)

        selC3 = u.select_atoms(f"name C3 and resid {straight_mol}").positions[0]
        selC14 = u.select_atoms(f"name C13 and resid {straight_mol}").positions[0]
        # sel_all = u.select_atoms("all")

        for_cross_prod = [selC14[0]-selC3[0], selC14[1]-selC3[1], selC14[2]-selC3[2]]
        for_cross_prod = for_cross_prod/np.linalg.norm(for_cross_prod)
        # for_cross_prod = pff.Avg_Vector(u)

        cross_prod = np.cross(for_cross_prod, Zvec)
        cross_prod_norm = cross_prod/(np.sqrt(cross_prod[0]**2 + cross_prod[1]**2 + cross_prod[2]**2))
        rotation_angle = pff.Vector_Plane_Angle(u, Zvec)[0]
        ag = u.atoms
        
        rotated = rotateby(rotation_angle*rot_dir, cross_prod_norm, ag=ag)(ag)
        rotated.write(f'{no_file_extensnion}_{axis_of_rotation}_rotated.gro', pbc=True)
        print(rotation_angle, cross_prod_norm)
    else:
        u = mda.Universe(input_top, input_trr)

        with mda.Writer(f'{no_file_extensnion}_{axis_of_rotation}_rotated.trr', n_atoms = u.atoms.n_atoms, multiframe = True) as W:
            
            try:
                if end == -2 or end == len(u.trajectory):
                    traj = u.trajectory[begin::stride]
                else:
                    traj = u.trajectory[begin : end + 1: stride]
            except IndexError:
                print("index error\nexiting...")
                exit(7)
            
            for ts in traj:#u.trajectory[:1000]:
                selC3 = u.select_atoms(f"name C3 and resid {straight_mol}").positions[0]
                selC14 = u.select_atoms(f"name C13 and resid {straight_mol}").positions[0]  

                for_cross_prod = [selC14[0]-selC3[0], selC14[1]-selC3[1], selC14[2]-selC3[2]]
                for_cross_prod_norm = for_cross_prod/np.linalg.norm(for_cross_prod) # normalize the molecule vector
                rotation_angle = np.rad2deg(np.arccos(np.dot(for_cross_prod_norm, Zvec))) 
                cross_prod = np.cross(for_cross_prod_norm, Zvec) # orthogonal vector to the molecule vector and Z axis vector

                # cross_prod_norm = cross_prod/(np.sqrt(cross_prod[0]**2 + cross_prod[1]**2 + cross_prod[2]**2))
                # rotation_angle = pff.Vector_Plane_Angle(u, Zvec)[0]
                ag = u.atoms

                rotated = rotateby(rotation_angle*rot_dir, cross_prod, ag=ag)(ag)
                u.atoms = rotated

                W.write(u.atoms)

if __name__ == "__main__":
    main(input_top= args.s, input_trr= args.f, begin = args.begin, end = args.end, stride = args.stride, straight_mol = args.straight_mol, axis_of_rotation = args.axis_of_rotation, rot_dir = args.rot_dir)


