# Performs comparison of Control and Stenosis generation resistances to determine if a segment contains severe to-be-fixed stenoses

import numpy as np
import json

from functions.utils import *
from functions.project_org import ProjectPath
from functions.find_generations import *
from functions.solver_io import Solver0D

########################
# Assistance Functions #
########################

def diameter_formula(order, age):
    ''' from https://journals.physiology.org/doi/full/10.1152/ajpheart.00123.2020 
    We assume order is roughly 16 - gen'''
    a = 1.203e-4
    b = 0.3927
    c = 0.2381
    d = 0.001 + a * order * np.exp(b * order ) * (age ** c)
    return d

def get_control_resistance(vessel_gen, model_age, vessel, viscosity = .04):

    # we assume vessel_order is approximately 16 - vessel_gen since gen 0 = mpa = order 16 
    vessel_order = 16 - vessel_gen
    rad = diameter_formula(vessel_order, model_age) /2
    length = vessel['vessel_length']
    control_resistance = 8 * viscosity * length / (np.pi * (rad ** 4))
    return control_resistance
    
    
def get_generation_resistance_averages(control_generations):
    ''' computes average of resistances over a particular generation by the total vessel length of that generation'''
    resistance_averages = {}
    for gen, branches in control_generations.items():
        total_resistance = 0
        total_length = 0
        for branch, vessel_segs in branches.items():
            for vessel in vessel_segs:
                total_length = vessel['vessel_length']
                total_resistance += vessel['vessel_length'] * vessel['zero_d_element_values']['R_poiseuille']
        
        resistance_averages[gen] = total_resistance / total_length

    return resistance_averages

def find_stenosis_vessels(test_gens, test_age,  r_threshold, plausible_gens = {0,1,2,3}):
    ''' finds stenosis vessels in test_gens
    r_threshold is the percent higher the test resistance is than the control,
    This should only hold true for the first 4 generations = orders 16 - 12'''
    total_len = 0
    sten_len = 0
    stenosis_vessels = []
    control_vessel_radii = []
    for gen, branches in test_gens.items():
        for branch, vessel_segs in branches.items():
            for vessel in vessel_segs:
                total_len += vessel['vessel_length']
                if gen in plausible_gens:
                    # if the r poiseuille values is greater than the average
                    control_res = get_control_resistance(gen, test_age, vessel, viscosity=.04 )
                    if vessel['zero_d_element_values']['R_poiseuille'] > ( r_threshold * control_res):
                        print(vessel['vessel_id'],vessel['zero_d_element_values']['R_poiseuille'],get_control_resistance(gen, test_age, vessel, viscosity=.04 ), vessel['zero_d_element_values']['stenosis_coefficient'])
                        control_vessel_radii.append(diameter_formula((16 - gen), test_age) / 2)
                        stenosis_vessels.append(vessel['vessel_id'])
                        sten_len += vessel['vessel_length']
    return stenosis_vessels, control_vessel_radii, total_len, sten_len
                    
        

########
# Main #
########

def tool_main(args):
    raise NotImplemented

def dev_main(args):
    
    org = ProjectPath(args.root)
    
    ## go through each specified model
    for model_name in args.models:
        test_model = org.find_model(model_name)
        if test_model.type != 'nci_stenosis':
            print('Model does not belong to nci_stenosis type.')
        else:
            
            ## retrieve test model
            test_age = int(test_model.params['metadata']['age'])
            
            # get test generations
            test_solver = Solver0D()
            test_solver.read_solver_file(test_model.model_solver)
            test_gens = find_generations(test_solver)
            
            
            r_threshold = args.r_threshold # 4x  higher
            plausible_gens  = {0,1,2,3} # first 4 generations
            vessel_ids, control_vessel_radii, total_len, sten_len = find_stenosis_vessels( test_gens, test_age, r_threshold, plausible_gens) 
            
            print(vessel_ids)
            with open(test_model.stenosis_vessels_file, 'w') as sten_file:
                json.dump({
                           'r_threshold': r_threshold ,
                           'stenosis_vessel_ids': vessel_ids,
                           'control_vessel_radii': control_vessel_radii,
                           'total_tree_len': total_len,
                           'total_sten_length': sten_len},
                          sten_file,  indent = 4, sort_keys=True)
        
            
            
            

if __name__ == '__main__':
    
    parser, _, dev, tool = create_parser(description='Generates artificial stenosis files')
    
    
    # dev params
    dev.add_argument('-root', dest = 'root', type = str, default = '.',  help = 'Root to entire project')
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-r_threshold', type = int, default = 4,  help = 'how many x control resistance the stenosis resistance must be to be considered a fixable location')
    args = parser.parse_args()
    
    
        
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)