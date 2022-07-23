from src.solver import Solver0D
from src.run_sim import run_sim, validate_rez, get_waveform_file
import argparse
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', dest = 'file', help = 'solver file to run')
    parser.add_argument('-s', dest = 'csv', action = 'store_true', default = False, help = 'save csv file: Default = False')
    parser.add_argument('-b', dest = 'branch', action = 'store_true', default=False,  help = 'to convert c output to python branch files: Default = False')
    parser.add_argument('-nv', dest = 'validate', action = 'store_false', default = True, help = 'validate the run with inlet pressure waveform: Default = True')
    parser.add_argument('-nss', dest = 'steady_sol', action = 'store_false', default = True, help = 'use steady solutions: Default = True')
    

    
    args = parser.parse_args()
    
    solver = Solver0D()
    solver.read_solver_file(args.file)
    rez = run_sim(solver = solver, use_steady_soltns=args.steady_sol, save_branch_results=args.branch, save_csv = args.csv, debug = True)
    if args.validate:
        print('Saving waveforms...', end = '\t', flush  =True)
        validate_rez(solver=solver, sim_results=rez, out_file = get_waveform_file(solver.solver_file))
        print('Done')
    
    
    