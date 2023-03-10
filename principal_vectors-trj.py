'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023, 5, 2327. doi: 10.3390/molecules28052327
'''
from importlib.metadata import metadata
import MDAnalysis as mda 
from MDAnalysis.transformations import rotateby
import numpy as np
import pandas as pd
import Plane_fitting_functions as pff
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from matplotlib.animation import FFMpegWriter
import argparse

parser = argparse.ArgumentParser(description='Imaging of molecules to reconstruct their real-space coordinates.')
parser.add_argument('-f', '--f', type = str, metavar = '', required = True, help = 'Trajectory file')
parser.add_argument('-s', '--s', type = str, metavar = '', required = True, help = 'Topology file')
parser.add_argument('-s2', '--s2', type = str, metavar = '', required = True, help = 'File with coordinates of the system where the molecules are oriented along the z axis with their long axis.')
parser.add_argument('-stm', '--straight_mol', type = str, metavar = '', required = True, help = 'Residue ID of a straight molecule (all trans configurations)')
parser.add_argument('-rd', '--rot_dir', type = int, metavar = '', choices=[1, -1], required = False, help = 'Clockwise/Anticlockwise rotation. 1 - Clockwise, -1 - Anticlockwise. Default is 1', 
        default = 1)
args = parser.parse_args()

def main(input_trr, input_top, input_top2, straight_mol, rot_dir):
    no_file_extensnion = str(input_trr[:len(input_trr) - 4])

    u = mda.Universe(input_top, input_trr)
    u2 = mda.Universe(input_top)  
    u3 = mda.Universe(input_top2)
    ''' Rotate the topology file to allign with the z axis '''
    selC3_last_frame = u2.select_atoms(f"name C3 and {straight_mol}").positions[0]
    sel_C13_last_frame = u2.select_atoms(f"name C13 and {straight_mol}").positions[0]

    for_cross_prod = [sel_C13_last_frame[0]-selC3_last_frame[0], sel_C13_last_frame[1]-selC3_last_frame[1], sel_C13_last_frame[2]-selC3_last_frame[2]]
    for_cross_prod = for_cross_prod/np.linalg.norm(for_cross_prod)

    Zvec = [0, 0, 1]
    r_angle = pff.Vector_Plane_Angle(u2, Zvec)[0]
    cross_prod = np.cross(for_cross_prod, Zvec)
    ''' done '''

    sel_C8 = u3.select_atoms("name C8").positions
    intercept_time_x = sel_C8.T[0]/10
    intercept_time_y = sel_C8.T[1]/10

    num_col = int(len(u.select_atoms('name C8').residues.resids))
    x_phi_traj = np.empty((0, num_col))
    y_phi_traj = np.empty((0, num_col))

    for ts in u.trajectory[ : : 20]:
        
        '''System Rotation'''
        ag = u.atoms
        rotated = rotateby(r_angle*rot_dir, cross_prod, ag=ag)(ag)
        '''System Rotation'''
        # rotated.atoms.write(f'{no_file_extensnion}_frame{ts.time}.gro')
        
        residues = u.select_atoms('name C8').residues.resids

        x_phi_I2 = []
        y_phi_I2 = []
        for i in residues:
            res_sel = rotated.select_atoms(f'resid {i}')

            p1, p2, p3 = abs(res_sel.principal_axes())
            # p1, p2, p3 = res_sel.principal_axes()
            cos_phi = np.arccos(p2[0])
            cos_phi = np.round(cos_phi, 2)
            x_phi_I2.append(cos_phi)

            sin_phi = np.arcsin(p2[0])
            sin_phi = np.round(sin_phi, 2)
            y_phi_I2.append(sin_phi)
        
        x_phi_I2 = np.array(x_phi_I2).reshape(-1, num_col)
        y_phi_I2 = np.array(y_phi_I2).reshape(-1, num_col)
        
        x_phi_traj = np.append(x_phi_traj, x_phi_I2, axis=0)
        y_phi_traj = np.append(y_phi_traj, y_phi_I2, axis=0)

    cm = 1/2.54
    pixel = px = 1/plt.rcParams['figure.dpi']
    # fig, ax = plt.subplots(figsize=(1920*pixel, 1080*pixel), tight_layout=True)
    fig, ax = plt.subplots()
    ax.set_xlabel('X direction', fontweight='bold')
    ax.set_ylabel('Y direction', fontweight='bold')
    ax.set_title('Second principal axis', fontweight='bold', color='blue')
    plt.xticks([])
    plt.yticks([])
    #ax.spines['left', 'right', 'top', 'bottom'].set_linewidth(0.2) #['left', 'right', 'top', 'bottom']

    l, = plt.plot([], [], 'k-o') # no idea ||| label='Bulk HEX crystallite 2'
    plt.tight_layout()
    plt.legend(loc='best', frameon=False, handlelength=0, title='Rotator, layer 2, 289K')
    metadata = dict(title=f'{no_file_extensnion}_p2', artist = 'S.I.')
    writer = FFMpegWriter(fps=20, metadata=metadata)

    with writer.saving(fig, f'{no_file_extensnion}_p2.mp4', dpi=200):
        for i in range(len(x_phi_traj)):
            x_data = intercept_time_x
            y_data = intercept_time_y
            U = x_phi_traj[i]/np.sqrt((x_phi_traj[i])**2 + (y_phi_traj[i]**2))
            V = y_phi_traj[i]/np.sqrt((x_phi_traj[i])**2 + (y_phi_traj[i]**2))
            Q = plt.quiver(x_data, y_data, U, V, pivot='middle', width=0.0045)
            writer.grab_frame()
            Q.remove()

if __name__ == '__main__':
    main(input_top= args.s, input_top2=args.s2, input_trr= args.f, straight_mol = args.straight_mol, rot_dir = args.rot_dir)