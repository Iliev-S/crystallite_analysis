'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023, 5, 2327. doi: 10.3390/molecules28052327
'''
from cmath import nan
import MDAnalysis as mda 
import numpy as np
from typing import List, Tuple, Dict 
from sortedcontainers import SortedDict
import colorama
from time import time

class Covariance_Matrix():
    def __init__(self, coordinates) -> None:
        self.coord = coordinates
        self.X = self.coord.T[0]
        self.Y = self.coord.T[1]
        self.Z = self.coord.T[2]

    def cov_mat(self):
        mat = np.array([self.X, self.Y, self.Z])
        cov = np.cov(mat)
        return cov

class PCAPlane(Covariance_Matrix):
    def __init__(self, coordinates) -> None:
        super().__init__(coordinates)

    def com(self):
        com = np.array([np.average(self.X), np.average(self.Y), np.average(self.Z)])
        return com

    def normal_vector(self):
        eigval, eigvec = np.linalg.eigh(self.cov_mat())
        return eigvec.T[0]
    
    def plane(self):
       normal = self.normal_vector()
       d = np.dot(normal, self.com())
       plane = np.round(np.array([normal[0], normal[1], normal[2], d]), 3)
       return plane


class PointOfIntercept:
    
    def __init__(self, p1 = None, p2 = None, plane = None) -> None:
        self.p1 = p1
        self.p2 = p2
        self.plane = plane

    def vector(self):
        return self.p2 - self.p1
    
    def vector_param(self):
        a, b, c, d =  self.plane[0], self.plane[1], self.plane[2], self.plane[3]
        px, py, pz = self.p1[0], self.p1[1], self.p1[2]
        vecx, vecy, vecz = self.vector()[0], self.vector()[1], self.vector()[2]  
        t = (d - a*px - b*py - c*pz)/(a*vecx + b*vecy + c*vecz)

        return t

    def intercept(self):
        px, py, pz = self.p1[0], self.p1[1], self.p1[2]
        vecx, vecy, vecz = self.vector()[0], self.vector()[1], self.vector()[2]  
        t = self.vector_param()
        x = px + t*vecx
        y = py + t*vecy
        z = pz + t*vecz
        point_of_intercept = np.array([x, y, z])
        
        return point_of_intercept

    def rot_matrix(self):
        a, b, c, d =  self.plane[0], self.plane[1], self.plane[2], self.plane[3]
        plane_normal = np.array([a, b, c])
        xy_normal = np.array([0, 0, 1])
        try:
            c = np.dot(plane_normal, xy_normal)/(np.linalg.norm(plane_normal)*np.linalg.norm(xy_normal)) #costheta
            s = np.sqrt(1-c*c) #sintheta
            C = 1 - c
            direction = np.cross(plane_normal, xy_normal)/np.linalg.norm(np.cross(plane_normal, xy_normal))
            if np.isnan(direction).sum():
                print('initial point returned!')
                return self.intercept()
            x, y, z = direction[0], direction[1], direction[2]
            rmat_z  = np.array([[x*x*C + c, x*y*C - z*s, x*z*C + y*s], 
                                [y*x*C + z*s, y*y*C + c, y*z*C - x*s],
                                [z*x*C - y*s, z*y*C + x*s, z*z*C + c]])
            rotated_point =  np.dot(rmat_z, self.intercept().T) 
            return rotated_point
        except ZeroDivisionError or RuntimeWarning:
            return self.intercept()

    def mat_of_rot(self):
        a, b, c, d =  self.plane[0], self.plane[1], self.plane[2], self.plane[3]
        plane_normal = np.array([a, b, c])
        xy_normal = np.array([0, 0, 1])
        c = np.dot(plane_normal, xy_normal)/(np.linalg.norm(plane_normal)*np.linalg.norm(xy_normal)) #costheta
        s = np.sqrt(1-c*c) #sintheta
        C = 1 - c
        direction = np.cross(plane_normal, xy_normal)/np.linalg.norm(np.cross(plane_normal, xy_normal))
        if np.isnan(direction).sum():
            print('Given plane is parallel to the xy plane!')
            return np.identity(3)
        x, y, z = direction[0], direction[1], direction[2]
        rmat_z  = np.array([[x*x*C + c, x*y*C - z*s, x*z*C + y*s], 
                            [y*x*C + z*s, y*y*C + c, y*z*C - x*s],
                            [z*x*C - y*s, z*y*C + x*s, z*z*C + c]])
        return rmat_z        

class Vacuum_Layer:
    def __init__(self, positions, resids, criterion):
        self.positions = positions
        self.resids = resids
        self.criterion = criterion
        
    def hash_table_pos_res(self, ordered = True):
        pos_ind_hash: Dict[np.float32, int] = dict()
        for i, j in zip(self.positions, self.resids):
            pos_ind_hash[i] = int(j)
        if ordered is True:
            return SortedDict(pos_ind_hash)
        else:
            return pos_ind_hash
        
    def vacuum(self):
        empty_list: List[int] = list()
        empty_dict: Dict[np.float32, List[int]] = dict()
    
        s = self.hash_table_pos_res()
        start = 0
        for key in range(len(list(s.keys())) - 1):
            key2 = key + 1
            if s.keys()[key2] - s.keys()[key] > self.criterion:
                for i in range(start, key2):
                    empty_list.append(s.values()[i])
                start = i + 1
                empty_dict[s.keys()[key]] = empty_list
                empty_list = []
        if start - 1 != len(list(s.keys())):
            for j in range(start, key2+1):
                empty_list.append(s.values()[j])
            empty_dict[s.keys()[start]] = empty_list
            empty_list = []    
        return empty_dict

def Vector_Plane_Angle(universe: mda.Universe, plane_vector_normal: List[float]) -> Tuple[float, float, List[float]]:
    """
        Returns tuple of average, std, all_angles
    """

    # x, y, z coordinates of the C3 and C14 atoms in Hexadecane
    selection_C3 = universe.select_atoms("name C3").positions
    selection_C13 = universe.select_atoms("name C13").positions

    # list of C3-C13 vector coordinates
    C3C13 = [i-j for i,j in zip(selection_C13, selection_C3)]

    C3C13 = np.array(C3C13).flatten()

    # length of the plane normal
    plane_vec_length = np.sqrt(plane_vector_normal[0]**2 + plane_vector_normal[1]**2 + plane_vector_normal[2]**2)

    vec_plane_ang_rad = []

    for i in range(0, len(C3C13), 3):
        # scalar product between each vector in the given universe and the plane normal
        dotprod = abs(C3C13[i]*plane_vector_normal[0] + C3C13[i+1]*plane_vector_normal[1] + C3C13[i+2]*plane_vector_normal[2])
        # lenght of each C3-C14 vector in the given universe
        mol_vec_length = np.sqrt(C3C13[i]*C3C13[i] + C3C13[i+1]*C3C13[i+1] + C3C13[i+2]*C3C13[i+2])
        # angle in degrees between each C3-C14 vector and the plane normal
        vec_plane_ang_rad.append(dotprod/(mol_vec_length*plane_vec_length))

    for i, _ in enumerate(vec_plane_ang_rad):
        if vec_plane_ang_rad[i] <= -1:
            vec_plane_ang_rad[i] += 1e-6
        
        if vec_plane_ang_rad[i] >= 1:
            vec_plane_ang_rad[i] -= 1e-6
    import warnings
    warnings.filterwarnings("error")
    try: 
        vec_plane_ang_deg = np.arccos(vec_plane_ang_rad)*57.2958
    except RuntimeWarning:
        print(vec_plane_ang_rad)
    vec_plane_ang_deg_avg = sum(vec_plane_ang_deg)/len(vec_plane_ang_deg)
    vec_plane_ang_deg_std = np.std(vec_plane_ang_deg)
    
    return vec_plane_ang_deg_avg, vec_plane_ang_deg_std, vec_plane_ang_deg


def Head_Tail(universe: mda.Universe, atom1: str = 'C1', atom2: str = 'C16', avg_std = True, raw_data = True) -> Tuple[float, float, List[float]]:
    """
    Returns tuple of average, std, all_Head_Tail_distances

    avg_std = True: returns the average and standard deviation of the Head-Tail distances 
    raw_data = True: returns all distances in a list
    """

    # x, y, z coordinates of the C1 and C16 atoms in Hexadecane
    selection_C1 = universe.select_atoms(f"name {atom1}").positions
    selection_C16 = universe.select_atoms(f"name {atom2}").positions 

    # list of C1-C16 vector coordinates
    C1C16 = [i-j for i,j in zip(selection_C16, selection_C1)]

    HT_dist = []

    for i in range(len(C1C16)):
        # lenght of each C1-C16 vector in the given universe
        HT_dist.append(np.sqrt(C1C16[i][0]*C1C16[i][0] + C1C16[i][1]*C1C16[i][1] + C1C16[i][2]*C1C16[i][2]))

    if avg_std is True and raw_data is True:
        HT_dist_av = sum(HT_dist)/len(HT_dist)
        HT_dist_std = np.std(HT_dist)
        return HT_dist_av, HT_dist_std, HT_dist 

    if avg_std is False and raw_data is True:
        return HT_dist 
    
    if avg_std is True and raw_data is False:
        HT_dist_av = sum(HT_dist)/len(HT_dist)
        HT_dist_std = np.std(HT_dist)
        return HT_dist_av, HT_dist_std 
        

def Avg_Vector(universe: mda.Universe, atom1: str = 'C3', atom2: str = 'C14') -> List[float]: 
    """ Returns a list with x, y, z coordinates of an averged normalized vector"""
    selection1 = universe.select_atoms(f'name {str(atom1)}').positions
    selection2 = universe.select_atoms(f'name {str(atom2)}').positions

    sel1sel2 = [i-j for i,j in zip(selection1, selection2)]

    vectors = np.array(sel1sel2)
    X = vectors.T[0].flatten()
    Y = vectors.T[1].flatten()
    Z = vectors.T[2].flatten()

    avg_vector = [np.sum(X)/len(X), np.sum(Y)/len(Y), np.sum(Z)/len(Z)]
    avg_vector_norm = avg_vector/np.linalg.norm(avg_vector)

    return avg_vector_norm

def res_ind_dict(u) -> Dict:
    ''' A dictionary that maps each molecule in the system to its starting index'''
    residues = u.atoms.residues.resids
    index_list = []
    x = 0
    for i in range(len(residues)):
        index_list.append(x)
        x += len(u.select_atoms(f'resid {residues[i]}'))
    dictionary = {key:value for (key,value) in zip(residues, index_list)}
    return dictionary

def most_values(dictionary, num_of_layers = 1):
    def sorting_key(item):
        return len(item[1])
    if num_of_layers > 0  and num_of_layers <= len(dictionary.keys()):
        new_dict = dict(sorted(dictionary.items(), key=sorting_key)[::-1])
        return list(new_dict.keys())[:num_of_layers]
    else: 
        raise TypeError('Invalid number of keys!')

class Progression_timer():
    def __init__(self, progress, total) -> None:
        self.progress = progress
        self.total = total

    def run_once(f):
        def wrapper(*args, **kwargs):
            if not wrapper.has_run:
                wrapper.has_run = True
                return f(*args, **kwargs)
        wrapper.has_run = False
        return wrapper
    
    global __start_time
    __start_time = time()

    def time_estimator(self,__start_time, step=1):
        if self.progress%step == 0:
            elapsed_time = time() - __start_time
            estimated_time = elapsed_time * self.total / (self.progress+1) / 60
            #print(colour + f'{estimated_time:.2f} min')
            return estimated_time

    def progress_bar(self, colour=colorama.Fore.LIGHTBLACK_EX):
        percent = 100*(self.progress/float(self.total))
        timer = self.time_estimator(__start_time)
        remaining = timer*(100 - percent)/100
        bar = '█'*int(percent) + '-' * (100 - int(percent))
        if timer != None:
            print(colour + f'\r|{bar}| {percent:.2f}% | Rem. Time: {remaining:.2f} min, Tot. Time: {timer:.2f} min', end='\r')
        print(colour + f'\r|{bar}| {percent:.2f}%', end='\r')
        if self.progress == self.total:
            print(colorama.Fore.GREEN + f'\r|{bar}| {percent:.2f}% | Rem. Time: {remaining:.2f} min, Tot. Time: {timer:.2f} min', end='\r')
            print(colorama.Fore.RESET)

def layer_split(increment, start, end):
    result_list = []

    for j in range(increment):
        empty_list = []
        for i in range(start + j, end+1, increment):
            empty_list.append(i)
        result_list.append(empty_list)
    return result_list

# def most_values(dictionary):
#     maxcount = max(len(v) for v in dictionary.values())
#     return[k for k,v in dictionary.items() if len(v) == maxcount]

