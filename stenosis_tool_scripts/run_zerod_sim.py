from functions.run_sim import *
    
if __name__ == '__main__':
    
    
    run_sim(solver_file = sys.argv[1], branch = True, block_print = False)
    validate_rez(solver_file=sys.argv[1], waveform_name=f'{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_waveforms.png')