from functions.comparisons_org import ComparisonPath
from functions.solver_io import Solver0D

import os

if __name__ == '__main__':
    
    ''' For Comparing a new technique of stenosis junction handling '''
    
    root = '.'
    
    comporg = ComparisonPath(root)
    
    for model in comporg.models:
        model_name = model.params['metadata']['name']
        solver_file=model.model_solver
        solver = Solver0D()
        solver.read_solver_file(solver_file = solver_file)
        
        
        # for each junction, compute a stenosis coefficient to add to the outlet branches downstream
        for junc in solver.junctions:
            if junc['junction_type'] != 'internal_junction':
                print(junc['junction_name'])
                
                inlet = junc['inlet_vessels']
                outlets = junc['outlet_vessels']
                S0 = junc['areas'][0]
                s_outlets = junc['areas'][1:]
                
                assert len(outlets) == len(s_outlets), 'Number of areas does not match number of outlets'
                
                outlet_areas = list(zip(outlets, s_outlets))
                
                density = solver.simulation_params['density']
                for out, S1 in outlet_areas:
                    out_vess = solver.get_vessel(out)
                    
                    # update with stenosis coefficient of junction
                    out_vess['stenosis_coefficient_junction'] = 1.52 * density * ( (S0/S1 - 1) **2) / (2 * S0**2)
                    out_vess['zero_d_element_values']['stenosis_coefficient'] += out_vess['stenosis_coefficient_junction']
                    
            new_solver_file = os.path.join(os.path.join(comporg.dir_names[model_name], comporg.NEW_0D), model_name + '_new_model.in' )
            solver.write_solver_file(new_solver_file)
                
        
        