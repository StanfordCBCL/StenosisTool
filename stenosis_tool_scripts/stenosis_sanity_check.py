# checks if a stenosis is fixed, what would happen


import sys
import shutil
from functions.run_sim import *
from functions.solver_io import *
import json
import numpy as np
import os

def compute_radii(viscosity, length, rp):
    r = ((8 * viscosity * length) / (rp * np.pi)) ** (1/4)
    return r


def new_inductance(old_ind, rad_rat):
    
    return old_ind / (rad_rat ** 2)

def new_capacitance(old_c, rad_rat):
    return rad_rat**2 * old_c

def new_sten_coeff(old_sten_coeff, rad_rat):
    
    # modifications to a_0 and a_s
    return old_sten_coeff / (rad_rat**4)
    

def new_r_poiseuille(old_r_p, rad_rat):
    return old_r_p / (rad_rat ** 4)
    
    



if __name__ == '__main__': 
    solver_file = sys.argv[1]
    stenosis_file = sys.argv[2]
    outdir  = sys.argv[3]
    
    solver = Solver0D()
    solver.read_solver_file(solver_file)
    
    with open(stenosis_file, 'r') as sfile:
        sten_info = json.load(sfile)
    
    vess = sten_info['stenosis_vessel_ids']
    normal_radii = sten_info['control_vessel_radii']
    
    total_len = 0
    for vidx in range(len(vess)):
        vid = vess[vidx]
        norm_rad = normal_radii[vidx]
        avessel = solver.get_vessel(id = vid)
        sten_rad = compute_radii(solver.simulation_params['viscosity'], avessel['vessel_length'], avessel['zero_d_element_values']['R_poiseuille'])
        rad_rat = norm_rad/sten_rad
        
        avessel['zero_d_element_values']['L'] = new_inductance(avessel['zero_d_element_values']['L'], rad_rat) / 10
        avessel['zero_d_element_values']['C'] = new_capacitance(avessel['zero_d_element_values']['C'], rad_rat) / 10
        avessel['zero_d_element_values']['R_poiseuille'] = new_r_poiseuille(avessel['zero_d_element_values']['R_poiseuille'], rad_rat) / 10
        avessel['zero_d_element_values']['stenosis_coefficient'] = new_sten_coeff(avessel['zero_d_element_values']['stenosis_coefficient'], rad_rat)
        total_len += avessel['vessel_length']

    print(total_len)
    
    out_filename = f"{solver.simulation_params['model_name']}_fixed_stenosis.in"
    out_filepath = os.path.join(outdir, out_filename)
    solver.write_solver_file(out_filepath)
    
    
    #run a simulation
    run_sim(out_filepath)
    validate_rez(out_filepath, os.path.splitext(out_filename)[0] + '_pressure_waveforms.png')
        
