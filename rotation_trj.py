import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np 
from MDAnalysis.transformations import rotateby
import MDAnalysis.transformations as trans
import Plane_fitting_functions as pff


input_top = str('R2-2Serfices-1K-400-500ns_whole_cryst6.gro')
input_trr = str('')
no_file_extensnion = str(input_top[:len(input_top) - 4])

if input_trr == 'None' or input_trr == '':
    u = mda.Universe(input_top)

    selC3 = u.select_atoms("name C3 and resid 523").positions[0]
    selC14 = u.select_atoms("name C13 and resid 523").positions[0]
    # sel_all = u.select_atoms("all")

    for_cross_prod = [selC14[0]-selC3[0], selC14[1]-selC3[1], selC14[2]-selC3[2]]
    for_cross_prod = for_cross_prod/np.linalg.norm(for_cross_prod)
    # for_cross_prod = pff.Avg_Vector(u)

    Zvec = [0, 0, 1]
    cross_prod = np.cross(for_cross_prod, Zvec)
    cross_prod_norm = cross_prod/(np.sqrt(cross_prod[0]**2 + cross_prod[1]**2 + cross_prod[2]**2))
    rotation_angle = pff.Vector_Plane_Angle(u, Zvec)[0]
    ag = u.atoms
    
    rotated = rotateby(rotation_angle, cross_prod_norm, ag=ag)(ag)
    rotated.write(no_file_extensnion + '_Z_rotated.gro', pbc=True)
    print(rotation_angle, cross_prod_norm)
else:
    u = mda.Universe(input_top, input_trr)

    with mda.Writer(no_file_extensnion + "_Z_rotated.trr", n_atoms = u.atoms.n_atoms, multiframe = True) as W:
        
        for ts in u.trajectory[:1000]:
            selC3 = u.select_atoms("name C3 and resid 76").positions[0]
            selC14 = u.select_atoms("name C13 and resid 76").positions[0]  

            for_cross_prod = [selC14[0]-selC3[0], selC14[1]-selC3[1], selC14[2]-selC3[2]]
            for_cross_prod_norm = for_cross_prod/np.linalg.norm(for_cross_prod) # normalize the molecule vector
            Zvec = [1, 0, 0]
            rotation_angle = np.rad2deg(np.arccos(np.dot(for_cross_prod_norm, Zvec))) 
            cross_prod = np.cross(for_cross_prod_norm, Zvec) # orthogonal vector to the molecule vector and Z axis vector

            # cross_prod_norm = cross_prod/(np.sqrt(cross_prod[0]**2 + cross_prod[1]**2 + cross_prod[2]**2))
            # rotation_angle = pff.Vector_Plane_Angle(u, Zvec)[0]
            ag = u.atoms

            rotated = rotateby(rotation_angle*1, cross_prod, ag=ag)(ag)
            u.atoms = rotated

            W.write(u.atoms)

            # print(cross_prod_norm, rotation_angle)

