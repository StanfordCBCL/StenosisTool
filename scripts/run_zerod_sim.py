from src.run_sim import run_sim, validate_rez
import os
import argparse
import sys
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', dest = 'file', help = 'solver file to run')
    parser.add_argument('-nb', dest = 'branch', action = 'store_false', default=True,  help = 'save branch files: Default = True')
    parser.add_argument('-p', dest = 'print', action = 'store_true', default = False, help = 'to allow prints: Default = False')
    parser.add_argument('-nss', dest = 'steady_sol', action = 'store_false', default = True, help = 'use steady solutions: Default = True')
    parser.add_argument('-lc', dest = 'last_cycle', action = 'store_true', default = False, help = 'use steady solutions: Default = False')
    
    args = parser.parse_args()
    bp = not args.print
    run_sim(solver_file = args.file, save_branches = args.branch, block_print = args.print, use_steady_soltns=args.steady_sol, last_cycle=args.last_cycle)
    validate_rez(solver_file=args.file, waveform_name=f'{os.path.splitext(os.path.basename(args.file))[0]}_waveforms.png')
    