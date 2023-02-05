from lib2to3.pytree import convert
from .file_io import read_json
import numpy as np
from scipy.stats.qmc import Sobol, scale
from .lpn import Solver0D
import copy

# Functions for modifying parameters

def new_inductance(old_ind, rad_rat):
    return old_ind / (rad_rat ** 2)

def new_capacitance(old_c, rad_rat):
    return rad_rat**2 * old_c

def new_sten_coeff(old_sten_coeff, rad_rat):
    # modifications to a_0 and a_s
    return old_sten_coeff / (rad_rat**4)
    
def new_r_poiseuille(old_r_p, rad_rat):
    return old_r_p / (rad_rat ** 4)
    
def convert_all(vess, rad_rat):
    vess["zero_d_element_values"]["C"] = new_capacitance(vess["zero_d_element_values"]["C"], rad_rat)
    
    vess["zero_d_element_values"]["L"] = new_inductance(vess["zero_d_element_values"]["L"], rad_rat)
    vess["zero_d_element_values"]["R_poiseuille"] = new_r_poiseuille(vess["zero_d_element_values"]["R_poiseuille"], rad_rat)
    vess["zero_d_element_values"]["stenosis_coefficient"]= new_sten_coeff(vess["zero_d_element_values"]["stenosis_coefficient"], rad_rat)



def rp_to_radius(vessel, mu):
    
    return (8 * mu * vessel['vessel_length'] / (np.pi * vessel['zero_d_element_values']['R_poiseuille'])) ** (1/4)
        