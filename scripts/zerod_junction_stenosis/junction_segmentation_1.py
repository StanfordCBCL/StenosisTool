from src.data_org import DataPath, JunctionStenosisResults
from src.solver import Solver0D
from src.misc import create_parser

import os


def dev_main(args):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        model_results = JunctionStenosisResults(args.root, model)
        
        solver_file=model.model_solver
        solver = Solver0D()
        solver.read_solver_file(solver_file = solver_file)
        
        
        # for each junction, compute a stenosis coefficient to add to the outlet branches downstream
        for junc in solver.junctions:
            if junc['junction_type'] != 'internal_junction':
                
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
                    
            new_solver_file = os.path.join(os.path.join(model_results.method_dir_1, os.path.basename(solver_file).split('.')[0] + '_method_1.in' ))
            solver.write_solver_file(new_solver_file)
                
        
        

if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'Implementing a new method of stenosis junction handling')
    
    args = parser.parse_args()
    
    if args.mode == 'tool':
        pass
    if args.mode == 'dev':
        dev_main(args)
    
    