from src.solver import Solver0D
from src.run_sim import run_sim, validate_rez, get_waveform_file, get_branch_results_file
from src.misc import create_tool_parser, get_solver_name
import os
    
if __name__ == '__main__':

    parser = create_tool_parser(desc = 'Run simulations')

    parser.add_argument('-r', dest = 'recursive',action = 'store_true', default = False, help = 'Whether to run all files found in that tree (excluding tuning solvers), or just the specified dir')
    parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Force a simulation to run even if results already exist')
    parser.add_argument('-s', dest = 'csv', action = 'store_true', default = False, help = 'save csv file: Default = False')
    parser.add_argument('-b', dest = 'branch', action = 'store_true', default=False,  help = 'to convert c output to python branch files: Default = False')
    parser.add_argument('-nv', dest = 'validate', action = 'store_false', default = True, help = 'validate the run with inlet pressure waveform: Default = True')
    parser.add_argument('-nss', dest = 'steady_sol', action = 'store_false', default = True, help = 'use steady solutions: Default = True')
    

    
    args = parser.parse_args()
    for solver_dir in args.solver_dirs:
        
        solvers = []
        if args.recursive:
            for tup in os.walk(solver_dir):
                sdir = tup[0]
                sfile_name = get_solver_name(sdir)
                if sfile_name is not None:
                    sfile = os.path.join(sdir, sfile_name)
                    if sfile_name != 'tune_solver.in': # as long as its not a tuning solver
                        solvers.append(sfile)
        else:
            solvers.append(os.path.join(solver_dir, get_solver_name(solver_dir)))
        print(solvers)
        for solver_file in solvers:
            if not args.force and os.path.exists(get_branch_results_file(solver_file, cpp = True)):
                print(f'{solver_file} has already been ran. Skipping.')
                continue
            else:
                print(f'Running Sim for {solver_file}')
        
                solver = Solver0D()
                solver.read_solver_file(solver_file)
                rez = run_sim(solver = solver, use_steady_soltns=args.steady_sol, save_branch_results=args.branch, save_csv = args.csv, save_last = False, debug = True)
                if args.validate:
                    print('Saving waveforms...', end = '\t', flush  =True)
                    validate_rez(solver=solver, sim_results=rez, out_file = get_waveform_file(solver.solver_file))
                    print('Done')
                
            
            
        