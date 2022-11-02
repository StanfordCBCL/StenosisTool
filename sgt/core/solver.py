from re import S
import numpy as np
import pandas as pd  
from .lpn import LPN
from svzerodsolver.runnercpp import run_from_config

class Solver0Dcpp():
    ''' Interface to 0D C++ Solver
    '''
    
    def __init__(self, lpn: LPN, use_steady = True, last_cycle_only = True, mean_only = False, debug = False):
    
        self.lpn = lpn
        self.lpn.simulation_params['steady_initial'] = use_steady
        self.lpn.simulation_params['output_last_cycle_only'] = last_cycle_only
        self.lpn.simulation_params['output_mean_only'] = mean_only
        self.debug = debug

    def _print(self, s, end = '\n', flush = False):
        ''' debug print '''
        if self.debug:
            print(s, end = end, flush = flush)
    
    def run_sim(self):
        ''' run a simulation '''
        
        self._print("Running solver...", end = '\t', flush = True)
    
        results_df = run_from_config(self.lpn.lpn_data)
        self._print('Done')
        
        return SolverResults(results_df)
        

class SolverResults():
    ''' 0D C solver results file
    '''
    def __init__(self, df: pd.DataFrame):
        self.result_df = df
        
    def only_last_cycle(df, tc):
        ''' Returns a Solver Results with only last cycle
        '''
        df = df.copy()
        df = df[df['time'] >= df['time'].max() - tc].copy()
        df['time'] -= df['time'].min()
        return SolverResults(df)
       
    @classmethod
    def from_csv(cls, fp):
        ''' Load from a csv
        '''
        df = pd.read_csv(fp)
        return cls(df)
    
    def vessel_df(self, vessel_name):
        ''' retrieves a df isolated by name
        '''
        return self.result_df[self.result_df['name'] == vessel_name]
    
    def get_vessel_names(self):
        ''' get a list of all vessel names
        '''
        tmp = self.result_df.groupby('name').max()
        return list(tmp.index)
    
    def convert_to_python(self, lpn: LPN, out_file):
        ''' Convert the c results into python branch_results (not sorted)
        '''
        branch_results = {}
        
        # write time
        branch_results['time'] = np.array(list(self.vessel_df('V0')['time']))
        
        branch_results['distance'] = {}
        branch_results['flow'] = {}
        branch_results['pressure'] = {}
        
        branch_tree = lpn.get_branch_tree()

        for node in lpn.tree_bfs_iterator(branch_tree):
            branch_id = node.branch_id
            # get initial values for first vessel
            vess_Q = []
            vess_P = []
            vess_D = [0.0]
            vess_df = self.vessel_df('V' + str(node.vess_id[0]))
            vess_Q.append(list(vess_df['flow_in']))
            vess_Q.append(list(vess_df['flow_out']))
            vess_P.append(list(vess_df['pressure_in']))
            vess_P.append(list(vess_df['pressure_out']))
            vess_D.append(node.vessel_info[0]['vessel_length'])
            
            # if longer than just 1 vessel segment per branch
            if len(node.vess_id) > 1:
                cur_idx = 1
                for vess_id in node.vess_id[1:]:
                    vess_df = self.vessel_df('V' + str(vess_id))
                    vess_Q.append(list(vess_df['flow_out']))
                    vess_P.append(list(vess_df['pressure_out']))
                    vess_D.append(vess_D[-1] + node.vessel_info[cur_idx]['vessel_length'])
                    cur_idx += 1
        
            branch_results['distance'][branch_id] = np.array(vess_D)
            branch_results['flow'][branch_id] = np.array(vess_Q)
            branch_results['pressure'][branch_id] = np.array(vess_P)
        np.save(out_file, branch_results, allow_pickle=True)
    
    
    def save_csv(self, out_file):
        ''' save as a csv
        '''
        self.result_df.to_csv(out_file, sep = ',', header = True, index = False)
        

