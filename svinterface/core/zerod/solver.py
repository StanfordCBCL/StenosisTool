
import numpy as np
import pandas as pd  
from .lpn import LPN
from collections import defaultdict
from svzerodsolver.runnercpp import run_from_config
import matplotlib.pyplot as plt
from pathlib import Path
from svinterface.utils.misc import d2m
from svinterface.core.polydata import Centerlines
from scipy.interpolate import interp1d
from vtk.util.numpy_support import numpy_to_vtk as n2v

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

        # save csv
        if save_csv:
            print("Saving csv...", end = '\t', flush = True)
            results.save_csv(str(out_dir / 'branch_results.csv'))
            print("Done")

        # save branch results
        if save_branch:
            print("Converting to python branch results...", end = '\t', flush = True)
            branch_results = results.convert_to_python(self.lpn)
            np.save(str(out_dir / 'branch_results.npy'), branch_results, allow_pickle = True)
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
    
    def convert_to_mmHg(self):
        '''Performs conversion on all pressures to mmHg (will cause errors when applied multiple times)
        '''
        self.result_df['pressure_in'] = d2m(self.result_df['pressure_in'])
        self.result_df['pressure_out'] = d2m(self.result_df['pressure_out'])
    
    def validate_results(self, lpn: LPN, outfile, targets = None ):
        ''' plots the inlet pressure for last 3 cycles
        Assumes the entire solution was saved and not only last cycle
        
        targets if we want to include targets.
        '''
        # retrieve time of cardiac cycle and number of points per cycle
        inflow_tc = lpn.inflow.tc
        num_pts = int(lpn.simulation_params['number_of_time_pts_per_cardiac_cycle'])
        
        ## save flow and pressure graphs (last 3 cycles)
        fig, ax = plt.subplots(1,1 ,figsize=(15, 10))

        # plot last 3 cycles
        last_three_cycles = -3 *  num_pts
        
        # find mpa
        mpa = lpn.get_mpa_branch()
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
        
        # plot targets
        if targets:
            factor1 = 0
            if targets['maxPAP'][0] == targets['maxPAP'][1]:
                factor1 = 0.5
            ax.fill_between(x = time_ltc, y1 = targets['maxPAP'][0]-factor1, y2 = targets['maxPAP'][1]+factor1, color = 'r', alpha = .3, label = f"Target Max PAP: ({targets['maxPAP'][0]}, {targets['maxPAP'][1]})")
            factor2 = 0
            if targets['minPAP'][0] == targets['minPAP'][1]:
                factor2 = 0.5
            ax.fill_between(x = time_ltc, y1 = targets['minPAP'][0]-factor2, y2 = targets['minPAP'][1]+factor2, color = 'g', alpha = .3, label = f"Target Min PAP: ({targets['minPAP'][0]}, {targets['minPAP'][1]})")
            factor3 = 0
            if targets['mPAP'][0] == targets['mPAP'][1]:
                factor3 = 0.5
            ax.fill_between(x = time_ltc, y1 = targets['mPAP'][0]-factor3, y2 = targets['mPAP'][1]+factor3, color = 'b', alpha = .3, label = f"Target Avg PAP: ({targets['mPAP'][0]}, {targets['mPAP'][1]})")

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
    
    def convert_to_python(self, lpn: LPN):
        ''' Convert the c results into python branch_results (not sorted)
        '''
        branch_results = {}
        
        # write time
        branch_results['time'] = self.result_df['time'].unique()
        
        branch_results['distance'] = {}
        branch_results['flow'] = {}
        branch_results['pressure'] = {}
        
        branch_tree = lpn.get_tree()

        for node in lpn.tree_bfs_iterator(branch_tree, allow='branch'):
            branch_id = node.id
            # get initial values for first vessel
            vess_Q = []
            vess_P = []
            vess_D = [0.0]
            vess_df = self.vessel_df(node.vessel_info[0]['vessel_name'])
            vess_Q.append(list(vess_df['flow_in']))
            vess_Q.append(list(vess_df['flow_out']))
            vess_P.append(list(vess_df['pressure_in']))
            vess_P.append(list(vess_df['pressure_out']))
            vess_D.append(node.vessel_info[0]['vessel_length'])
            
            # if longer than just 1 vessel segment per branch
            if len(node.ids) > 1:
                cur_idx = 1
                for vessel in node.vessel_info[1:]:
                    vess_df = self.vessel_df(vessel['vessel_name'])
                    vess_Q.append(list(vess_df['flow_out']))
                    vess_P.append(list(vess_df['pressure_out']))
                    vess_D.append(vess_D[-1] + node.vessel_info[cur_idx]['vessel_length'])
                    cur_idx += 1
        
            branch_results['distance'][branch_id] = np.array(vess_D)
            branch_results['flow'][branch_id] = np.array(vess_Q)
            branch_results['pressure'][branch_id] = np.array(vess_P)
        return branch_results
    
    def save_csv(self, out_file):
        ''' save as a csv
        '''
        self.result_df.to_csv(out_file, sep = ',', header = True, index = False)
        
    def project_to_centerline(self, lpn: LPN, centerlines: Centerlines):
        """
        Project rom results onto the centerline
        Modified from: https://github.com/SimVascular/SimVascular/blob/master/Python/site-packages/sv_rom_extract_results/post.py
        """
        
        # convert c++ to python results
        results = self.convert_to_python(lpn)
        
        # assemble output dict
        rec_dd = lambda: defaultdict(rec_dd)
        arrays = rec_dd()

        # extract point arrays from geometries
        arrays_cent = {}
        for arr_name in centerlines.get_pointdata_arraynames():
            arrays_cent[arr_name] = centerlines.get_pointdata_array(arr_name)
        
        # add centerline arrays
        for name, data in arrays_cent.items():
            arrays[name] = data

        # centerline points
        points = centerlines.get_points()

        # all branch ids in centerline
        ids_cent = np.unique(arrays_cent['BranchId']).tolist()
        ids_cent.remove(-1)
        
        # retrieve time
        times = self.result_df['time'].unique()

        # loop all result fields
        for f in ['flow', 'pressure']:
            if f not in results:
                continue

            # check if ROM branch has same ids as centerline
            ids_rom = list(results[f].keys())
            ids_rom.sort()
            assert ids_cent == ids_rom, 'Centerline and ROM results have different branch ids'

            # initialize output arrays
            array_f = np.zeros((arrays_cent['Path'].shape[0], len(times)))
            n_outlet = np.zeros(arrays_cent['Path'].shape[0])

            # loop all branches
            for br in results[f].keys():
                # results of this branch
                res_br = results[f][br]

                # get centerline path
                path_cent = arrays_cent['Path'][arrays_cent['BranchId'] == br]

                # get node locations from 0D results
                path_0d_res = results['distance'][br]
                f_res = res_br
            
                assert np.isclose(path_0d_res[0], 0.0), 'ROM branch path does not start at 0'
                assert np.isclose(path_cent[0], 0.0), 'Centerline branch path does not start at 0'
                msg = 'ROM results and centerline have different branch path lengths'
                assert np.isclose(path_0d_res[-1], path_cent[-1]), msg

                # interpolate ROM onto centerline
                # limit to interval [0,1] to avoid extrapolation error interp1d due to slightly incompatible lenghts
                f_cent = interp1d(path_0d_res / path_0d_res[-1], f_res.T)(path_cent / path_cent[-1]).T

                # store results of this path
                array_f[arrays_cent['BranchId'] == br] = f_cent[:, range(len(times))]

                # add upstream part of branch within junction
                if br == 0:
                    continue

                # first point of branch
                ip = np.where(arrays_cent['BranchId'] == br)[0][0]

                # centerline that passes through branch (first occurence)
                cid = np.where(arrays_cent['CenterlineId'][ip])[0][0]

                # id of upstream junction
                jc = arrays_cent['BifurcationId'][ip - 1]

                # centerline within junction
                is_jc = arrays_cent['BifurcationId'] == jc
                jc_cent = np.where(np.logical_and(is_jc, arrays_cent['CenterlineId'][:, cid]))[0]
                
                # length of centerline within junction
                jc_path = np.append(0, np.cumsum(np.linalg.norm(np.diff(points[jc_cent], axis=0), axis=1)))
                jc_path /= jc_path[-1]
                
                # results at upstream branch
                res_br_u = results[f][arrays_cent['BranchId'][jc_cent[0] - 1]]

                # results at beginning and end of centerline within junction
                f0 = res_br_u[-1][range(len(times))]
                f1 = res_br[0][range(len(times))]
               
                # map 1d results to centerline using paths
                array_f[jc_cent] += interp1d([0, 1], np.vstack((f0, f1)).T)(jc_path).T

                # count number of outlets of this junction
                n_outlet[jc_cent] += 1

            # normalize results within junctions by number of junction outlets
            is_jc = n_outlet > 0
            array_f[is_jc] = (array_f[is_jc].T / n_outlet[is_jc]).T

            # assemble time steps
            for i, t in enumerate(times):
                arrays[f + '_' + str(t)] = array_f[:, i]
                
            # compute summary statistics
            avg = np.trapz( array_f, times, axis = 1) / (times[-1] - times[0])
            arrays['avg_' + f] = avg
            sys_tidx = np.argmin(abs(times - lpn.inflow.max_inflow_t))
            arrays['sys_' + f + '_' + str(times[sys_tidx])] = array_f[:, sys_tidx]
            dia_tidx = np.argmin(abs(times - lpn.inflow.min_inflow_t))
            arrays['dia_' + f + '_' + str(times[dia_tidx])] = array_f[:, dia_tidx]
            
        # add arrays to centerline and write to file
        for f, a in arrays.items():
            centerlines.add_pointdata(array=a, array_name=f)
        
        return centerlines