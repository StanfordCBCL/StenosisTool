
import numpy as np
import pandas as pd  
from .lpn import LPN
from svzerodsolver.runnercpp import run_from_config
import matplotlib.pyplot as plt
from pathlib import Path
from sgt.utils.misc import d2m

class Solver0Dcpp():
    ''' Interface to 0D C++ Solver
    '''
    
    def __init__(self, lpn: LPN, use_steady = True, last_cycle_only = True, mean_only = False, debug = False):
    
        self.lpn = lpn
        self.use_steady = use_steady
        self.last_cycle_only = last_cycle_only
        self.mean_only = mean_only
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
    
    def run_sim_pipeline(self, validate, save_csv, save_branch, out_dir):
            # resolve validation and last_cycle/mean_only conflict
        if validate and self.mean_only:
            print("Cannot validate with mean_only.")
            exit(1)
        elif validate and self.last_cycle_only:
            self.lpn.simulation_params['output_last_cycle_only'] = False
            post_convert_last_cycle = True
        else:
            post_convert_last_cycle = False

        out_dir = Path(out_dir)
    
        results = self.run_sim()

        # validate
        if validate:
            print("Validating results...", end = '\t', flush = True)
            results.validate_results(self.lpn, str(out_dir / "inlet_pressure_validation.png"))
            print("Done")
        
        # convert to last cycle
        if post_convert_last_cycle:
            results = results.only_last_cycle(self.lpn.inflow.tc)
            # reset back to avoid future issue
            self.lpn.simulation_params['output_last_cycle_only'] = True

    
        if save_csv:
            print("Saving csv...", end = '\t', flush = True)
            results.save_csv(str(out_dir / 'branch_results.csv'))
            print("Done")
    
        if save_branch:
            print("Converting to python branch results...", end = '\t', flush = True)
            results.convert_to_python(self.lpn, str(out_dir / 'branch_results.npy'))
            print("Done")
        return results

class SolverResults():
    ''' 0D C solver results file
    '''
    def __init__(self, df: pd.DataFrame):
        self.result_df = df
        
    def only_last_cycle(self, tc):
        ''' Returns a Solver Results with only last cycle
        '''
        df = self.result_df.copy()
        df = df[df['time'] >= df['time'].max() - tc].copy()
        df['time'] -= df['time'].min()
        return SolverResults(df)
    
    def validate_results(self, lpn: LPN, outfile ):
        ''' plots the inlet pressure for last 3 cycles
        Assumes the entire solution was saved and not only last cycle
        '''
        # retrieve time of cardiac cycle and number of points per cycle
        inflow_tc = lpn.inflow.tc
        num_pts = int(lpn.simulation_params['number_of_time_pts_per_cardiac_cycle'])
        
        ## save flow and pressure graphs (last 3 cycles)
        fig, ax = plt.subplots(1,1 ,figsize=(15, 10))

        # plot last 3 cycles
        last_three_cycles = -3 *  num_pts
        
        # find mpa
        mpa = lpn.get_mpa()
        mpa_name = mpa.vessel_info[0]['vessel_name']
    
        v0 = self.vessel_df(mpa_name)
        inlet_pressure = d2m(np.array(v0['pressure_in'][last_three_cycles:]))
        time_ltc = np.array(v0['time'][last_three_cycles:])
        # plot pressure curves
        ax.plot(time_ltc, inlet_pressure)
        ax.set_title(f"Inlet Pressure", fontdict={'fontsize': 24})
        ax.set_xlabel('time (s)', fontdict={'fontsize': 20})
        ax.set_ylabel('pressure (mmHg)', fontdict={'fontsize': 20})
        
        # get last cardiac cycle
        results = self.only_last_cycle(inflow_tc)
        v0_lc =results.vessel_df(mpa_name)
        # compute mean, max, min PAP
        stable_inlet_pressure = d2m(np.array(v0_lc['pressure_in']))
        time_lc = np.array(v0_lc['time'])
        mpap_sim = np.trapz(stable_inlet_pressure, time_lc) / inflow_tc
        
        x_max = time_lc[np.where(stable_inlet_pressure == stable_inlet_pressure.max())] + time_ltc[0] + 2 * inflow_tc

        x_min = time_lc[np.where(stable_inlet_pressure == stable_inlet_pressure.min())] + time_ltc[0] + 2 * inflow_tc
        
        # plot mean, max, min
        ax.plot(x_max, 
                stable_inlet_pressure.max(), 
                'r^', 
                label = 'Max PAP: ' + str(round(stable_inlet_pressure.max(),2)))
        ax.plot(x_min, 
                stable_inlet_pressure.min(), 
                'g^', 
                label = 'Min PAP: ' + str(round(stable_inlet_pressure.min(), 2)))
        ax.hlines(y = mpap_sim, xmin = time_ltc[0], xmax = time_ltc[-1], linewidth=1, color='b', label = 'Avg PAP: ' + str(round(mpap_sim, 2)) )
        
        ax.tick_params(axis="x", labelsize=16) 
        ax.tick_params(axis = 'y', labelsize=16)
        ax.legend(fontsize = 24, loc = 'upper left', framealpha = .5)

        fig.savefig(outfile)
    
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
    
    def get_avg_val(self,vessel_name, val = 'flow_in'):
        ''' get the average flow of a result '''
        assert val in {'flow_in', 'flow_out', 'pressure_in', 'pressure_out'}, "Must be one of flow_in, flow_out, pressure_in, pressure_out"
        df = self.vessel_df(vessel_name)
        flow = np.array(df[val])
        time = np.array(df['time'])
        return np.trapz(flow, time) / (time[-1] - time[0])

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
        

