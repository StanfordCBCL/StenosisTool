# File: compare_3D_0D.py
# File Created: Sunday, 31st July 2022 5:32:43 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 18th October 2022 9:14:12 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Computes error measurements for a 0D and a 3D model


from src.misc import create_tool_parser, get_solver_path, get_solver_name, get_basename
from src.file_io import check_exists, write_json
from src.polydata import Centerlines
from src.lpn import Solver0D
from src.solver_results import SolverResults
from src.misc import m2d, d2m

import os
import numpy as np

class ErrorTerms():
    
    def __init__(self, centerlines3d: Centerlines, solver0d: Solver0D, solver_results: SolverResults) -> None:
        
        # save solver/results
        self.centerlines3d = centerlines3d
        self.solver0d = solver0d
        self.solver_results = solver_results
        
        # retrieve data from centerlines
        self.t3d, self.PAP, self.Q = self._get_cent_data(centerlines3d)
        
        # convert to last cycle
        self.solver_results.result_df = SolverResults.only_last_cycle(self.solver_results.result_df, tc = self.solver0d.inflow.tc)
        self.result_df = self.solver_results.result_df
        
        # map 0d times onto 3D
        t0d = np.array(self.result_df[self.result_df['name'] == 'V0']['time'])
        self.t0d = np.array([t0d[abs(t0d - i).argmin()] for i in self.t3d])
        
        assert len(self.t3d) == len(self.t0d), '3D and 0D times do not match in length'
        


    @staticmethod
    def _get_cent_data(centerlines3d):
        # retrieve all data from centerlines
        arr = centerlines3d.get_pointdata_arraynames()
        Q = []
        PAP = []
        t3d = []
        for name in arr:
            if name.startswith('pressure_'):
                t3d.append(float(name.split('_')[-1]))
                PAP.append(centerlines3d.get_pointdata(name))
            if name.startswith('flow_'):
                Q.append(centerlines3d.get_pointdata(name))

        return np.array(t3d), np.array(PAP), np.array(Q)
        
    
    @staticmethod
    def get_vname(vess_id):
        return 'V' + str(vess_id)
        


class CapErrors(ErrorTerms):
    
    def __init__(self, centerlines3d: Centerlines, solver0d: Solver0D, solver_results: SolverResults) -> None:
        super().__init__(centerlines3d, solver0d, solver_results)
        
        # list of branch id's of caps & vessel ids
        self.ibranches, self.obranches, self.ivessels, self.ovessels = self._get_cap_branches_vessels()
        
        # list of node ids
        branch_ids = self.centerlines3d.get_pointdata(self.centerlines3d.PointDataFields.BRANCHID)
        self.inodes = [np.where(branch_ids == i)[0][0] for i in self.ibranches]
        self.onodes = [np.where(branch_ids == o)[0][-1] for o in self.obranches]
        
        # only 1 inlet
        assert len(self.ibranches) == 1 and len(self.ivessels) == 1 and len(self.inodes) == 1, 'More than 1 inlet detected'
        
        # process solver rez to only have relevant data
        self.Q0d = []
        self.PAP0d = []
        for t in self.t0d:
            tmp = self.result_df[self.result_df['time'] == t]
            tmp_arr_q = []
            tmp_arr_pap = []
            tmp_arr_pap.append(tmp.loc[tmp['name'] == self.get_vname(self.ivessels[0]), 'pressure_in'])
            tmp_arr_q.append(tmp.loc[tmp['name'] == self.get_vname(self.ivessels[0]), 'flow_in'])
            for v in self.ovessels:
                tmp_arr_pap.append(tmp.loc[tmp['name'] == self.get_vname(v), 'pressure_out'])
                tmp_arr_q.append(tmp.loc[tmp['name'] == self.get_vname(v), 'flow_out'])
            self.Q0d.append(np.array(tmp_arr_q).flatten())
            self.PAP0d.append(d2m(np.array(tmp_arr_pap).flatten())) # convert to mmHG
        
        self.Q0d = np.array(self.Q0d)
        self.PAP0d = np.array(self.PAP0d)

        self.sys_t, self.dia_t = self._get_sys_dia_t()
        
    def _get_sys_dia_t(self):
        sys_t = 0
        dia_t = 0
        arr_names = self.centerlines3d.get_pointdata_arraynames()
        for name in arr_names:
            if name.startswith('sys'):
                sys_t = float(name.split('_')[-1])
            elif name.startswith('dia'):
                dia_t = float(name.split('_')[-1])
        return sys_t, dia_t
                
        
    def _get_cap_branches_vessels(self):
        
        inlet_branches = []
        outlet_branches = []
        inlet_vessels = []
        outlet_vessels = []
        branch_tree = self.solver0d.get_branch_tree()
        # determine outlet/inlet branches (should be mutually exclusive)
        for node in self.solver0d.tree_bfs_iterator(branch_tree):
            if node.get_type() == 'outlet':
                outlet_branches.append(node.branch_id)
                outlet_vessels.append(node.vess_id[-1])
            elif node.get_type() == 'inlet':
                inlet_branches.append(node.branch_id)
                inlet_vessels.append(node.vess_id[-1])
        
        return inlet_branches, outlet_branches, inlet_vessels, outlet_vessels

    ###################
    # ERROR FUNCTIONS #
    ###################

    # Define error functions
    def ep_avg(self):
        
        all_nodes = np.array(self.inodes + self.onodes)
        p3d = np.zeros(len(all_nodes))
        p0d = np.zeros(len(self.ivessels + self.ovessels))

        for tidx in range(len(self.t3d)):
            
            cur_p3d = self.PAP[tidx][all_nodes]
            cur_p0d = self.PAP0d[tidx]
            np.add(p3d, cur_p3d, p3d)
            np.add(p0d, np.abs(cur_p0d - cur_p3d), p0d)
            

        return np.divide(p0d, p3d).mean()
            
            
    def ep_sys(self):
        
        all_nodes = np.array(self.inodes + self.onodes)
        p3d = np.zeros(len(all_nodes))
        p0d = np.zeros(len(self.ivessels + self.ovessels))
        
        for tidx in range(len(self.t3d)):  
            cur_p3d = self.PAP[tidx][all_nodes]
            np.add(p3d, cur_p3d, p3d)
        
        # sys tidx 
        sys_tidx = np.where(self.t3d == self.sys_t)[0][0]
        sys_p3d = self.PAP[sys_tidx][all_nodes]
        sys_p0d = self.PAP0d[sys_tidx]
        np.add(p0d, np.abs(sys_p0d - sys_p3d), p0d)
        
        return np.divide(p0d, p3d).mean() * len(self.t3d)
    
    def ep_dia(self):
        
        all_nodes = np.array(self.inodes + self.onodes)
        p3d = np.zeros(len(all_nodes))
        p0d = np.zeros(len(self.ivessels + self.ovessels))
        
        for tidx in range(len(self.t3d)):  
            cur_p3d = self.PAP[tidx][all_nodes]
            np.add(p3d, cur_p3d, p3d)
        
        # dia tidx 
        dia_tidx = np.where(self.t3d == self.dia_t)[0][0]
        dia_p3d = self.PAP[dia_tidx][all_nodes]
        dia_p0d = self.PAP0d[dia_tidx]
        np.add(p0d, np.abs(dia_p0d - dia_p3d), p0d)
        
        return np.divide(p0d, p3d).mean() * len(self.t3d)
        
        
    def eq_avg(self):
        
        all_nodes = np.array(self.onodes)
        q3d_max = self.Q[0][all_nodes]
        q3d_min = self.Q[0][all_nodes]
        q0d = np.zeros(len(self.ovessels))
        for tidx in range(len(self.t3d)):
            
            cur_q3d = self.Q[tidx][all_nodes]
            cur_q0d = self.Q0d[tidx][1:] # ignore inlet
            q3d_max = np.where(q3d_max < cur_q3d, cur_q3d, q3d_max)
            q3d_min = np.where(q3d_min > cur_q3d, cur_q3d, q3d_min)
            np.add(q0d, np.abs(cur_q0d - cur_q3d), q0d)
        
        return np.divide(q0d, np.subtract(q3d_max, q3d_min)).mean() / len(self.t3d)
        
    
            
    def eq_sys(self):
        
        all_nodes = np.array(self.onodes)
        q3d_max = self.Q[0][all_nodes]
        q3d_min = self.Q[0][all_nodes]
        q0d = np.zeros(len(self.ovessels))
        for tidx in range(len(self.t3d)):
            
            cur_q3d = self.Q[tidx][all_nodes]
            q3d_max = np.where(q3d_max < cur_q3d, cur_q3d, q3d_max)
            q3d_min = np.where(q3d_min > cur_q3d, cur_q3d, q3d_min)
        
        # sys tidx 
        sys_tidx = np.where(self.t3d == self.sys_t)[0][0]
        sys_q3d = self.Q[sys_tidx][all_nodes]
        sys_q0d = self.Q0d[sys_tidx][1:]
        np.add(q0d, np.abs(sys_q0d - sys_q3d), q0d)
        
        return np.divide(q0d, np.subtract(q3d_max, q3d_min)).mean()
    
    def eq_dia(self):
        
        all_nodes = np.array(self.onodes)
        q3d_max = self.Q[0][all_nodes]
        q3d_min = self.Q[0][all_nodes]
        q0d = np.zeros(len(self.ovessels))
        for tidx in range(len(self.t3d)):
            
            cur_q3d = self.Q[tidx][all_nodes]
            q3d_max = np.where(q3d_max < cur_q3d, cur_q3d, q3d_max)
            q3d_min = np.where(q3d_min > cur_q3d, cur_q3d, q3d_min)
        
        # dia tidx 
        dia_tidx = np.where(self.t3d == self.dia_t)[0][0]
        dia_q3d = self.Q[dia_tidx][all_nodes]
        dia_q0d = self.Q0d[dia_tidx][1:]
        np.add(q0d, np.abs(dia_q0d - dia_q3d), q0d)
        
        return np.divide(q0d, np.subtract(q3d_max, q3d_min)).mean()
    
class CapErrorsNew(ErrorTerms):
    
    def __init__(self, centerlines3d: Centerlines, solver0d: Solver0D, solver_results: SolverResults) -> None:
        super().__init__(centerlines3d, solver0d, solver_results)
        # list of branch id's of caps & vessel ids
        self.ibranches, self.obranches, self.ivessels, self.ovessels = self._get_cap_branches_vessels()
        
        # list of node ids
        branch_ids = self.centerlines3d.get_pointdata(self.centerlines3d.PointDataFields.BRANCHID)
        self.inodes = [np.where(branch_ids == i)[0][0] for i in self.ibranches]
        self.onodes = [np.where(branch_ids == o)[0][-1] for o in self.obranches]
        
        # only 1 inlet
        assert len(self.ibranches) == 1 and len(self.ivessels) == 1 and len(self.inodes) == 1, 'More than 1 inlet detected'
        
    def _get_cap_branches_vessels(self):
        
        inlet_branches = []
        outlet_branches = []
        inlet_vessels = []
        outlet_vessels = []
        branch_tree = self.solver0d.get_branch_tree()
        # determine outlet/inlet branches (should be mutually exclusive)
        for node in self.solver0d.tree_bfs_iterator(branch_tree):
            if node.get_type() == 'outlet':
                outlet_branches.append(node.branch_id)
                outlet_vessels.append(node.vess_id[-1])
            elif node.get_type() == 'inlet':
                inlet_branches.append(node.branch_id)
                inlet_vessels.append(node.vess_id[-1])
        
        return inlet_branches, outlet_branches, inlet_vessels, outlet_vessels

    @staticmethod
    def get_avg(t, y):
        return np.trapz(y, t) / (t[-1] - t[0])
        
    def ep_avg(self):
        
        # MSE of caps,   (avg0d - avg3d / avg3d) = 
        # avg3d per cap
        nodes = np.array(self.inodes + self.onodes)
        avg3d = np.trapz( self.PAP[:, nodes],self.t3d, axis = 0) / (self.t3d[-1] - self.t3d[0])

        print(avg3d)
        
        time = np.array(self.solver_results.vessel_df('V0')['time'])
        avg0d = [self.get_avg(time, np.array(self.result_df[self.result_df['name'] == 'V' + str(v)]['pressure_in'])) for v in self.ivessels] + [self.get_avg(time, np.array(self.result_df[self.result_df['name'] == 'V' + str(v)]['pressure_out'])) for v in self.ovessels]
        avg0d = d2m(np.array(avg0d)) # convert to numpy and mmHg
        
        print(avg0d)
        
        return np.divide(np.abs(avg0d - avg3d), avg3d).mean()
        
        
    def ep_sys(self):
        pass
    def ep_dia(self):
        pass
    def eq_avg(self):
        # MSE of caps,   (avg0d - avg3d / avg3d) = 
        # avg3d per cap
        nodes = np.array( self.onodes)
        avg3d = np.trapz( self.Q[:, nodes],self.t3d, axis = 0) / (self.t3d[-1] - self.t3d[0])
        
        time = np.array(self.solver_results.vessel_df('V0')['time'])
        avg0d = [self.get_avg(time, np.array(self.result_df[self.result_df['name'] == 'V' + str(v)]['flow_out'])) for v in self.ovessels]
        avg0d = np.array(avg0d) # convert to numpy and mmHg
        
        print(avg0d)
        
        return np.divide(np.abs(avg0d - avg3d), avg3d).mean()
        
        
    def eq_sys(self):
        pass
    def eq_dia(self):
        pass

        
class PDErrors(CapErrors):
    
    def __init__(self, centerlines3d: Centerlines, solver0d: Solver0D, solver_results: SolverResults) -> None:
        super().__init__(centerlines3d, solver0d, solver_results)

    def e_pd_avg(self):
        
        inodes = np.array(self.inodes )
        onodes = np.array(self.onodes)
        dp3d = np.zeros(len(self.onodes))
        dp0d = np.zeros(len(self.ovessels))

        for tidx in range(len(self.t3d)):
            
            cur_dp3d = np.abs(self.PAP[tidx][inodes] - self.PAP[tidx][onodes])
            cur_dp0d = np.abs(self.PAP0d[tidx][0] - self.PAP0d[tidx][1:])
            np.add(dp3d, cur_dp3d, dp3d)
            np.add(dp0d, np.abs(cur_dp0d - cur_dp3d), dp0d)
            

        return np.divide(dp0d, dp3d).mean()
    
    def e_pd_sys(self):
        
        inodes = np.array(self.inodes)
        onodes = np.array(self.onodes)
        p3d = np.zeros(len(onodes))
        p0d = np.zeros(len(self.ovessels))
        
        for tidx in range(len(self.t3d)):  
            cur_p3d = np.abs(self.PAP[tidx][inodes] - self.PAP[tidx][onodes])
            np.add(p3d, cur_p3d, p3d)
        
        # sys tidx 
        sys_tidx = np.where(self.t3d == self.sys_t)[0][0]
        sys_p3d = np.abs(self.PAP[sys_tidx][inodes] - self.PAP[sys_tidx][onodes])
        sys_p0d = np.abs(self.PAP0d[sys_tidx][0] -  self.PAP0d[sys_tidx][1:])
        np.add(p0d, np.abs(sys_p0d - sys_p3d), p0d)
        
        return np.divide(p0d, p3d).mean() * len(self.t3d)
    
    def e_pd_dia(self):
        
        inodes = np.array(self.inodes)
        onodes = np.array(self.onodes)
        p3d = np.zeros(len(onodes))
        p0d = np.zeros(len(self.ovessels))
        
        for tidx in range(len(self.t3d)):  
            cur_p3d = np.abs(self.PAP[tidx][inodes] - self.PAP[tidx][onodes])
            np.add(p3d, cur_p3d, p3d)
        
        # dia tidx 
        dia_tidx = np.where(self.t3d == self.dia_t)[0][0]
        dia_p3d = np.abs(self.PAP[dia_tidx][inodes] - self.PAP[dia_tidx][onodes])
        dia_p0d = np.abs(self.PAP0d[dia_tidx][0] -  self.PAP0d[dia_tidx][1:])
        np.add(p0d, np.abs(dia_p0d - dia_p3d), p0d)
        
        return np.divide(p0d, p3d).mean() * len(self.t3d)
    
    
    
    
# Main
def main(args):
    
    for solver_dir in args.solver_dirs:
        print(f'Computing errors for {solver_dir}...')
        ## get 0D solver file
        solver_file = get_solver_path(solver_dir)
        solver0d = Solver0D()
        solver0d.read_solver_file(solver_file)
        
        ## get 0D solver results
        solver_rez_file = check_exists(os.path.join(solver_dir, get_basename(solver_file)+ '_branch_results.csv'), err = f"{os.path.join(solver_dir, get_basename(solver_file)+ '_branch_results.csv')} could not be found.", mkdir = False, ignore = False)
        solver0d_results = SolverResults.load_from_csv(solver_rez_file)
        
        ## load 3D file
        solver3d_results = Centerlines()
        three_d_dir = os.path.join(solver_dir, 'three_d_dir')
        solver3d_file = check_exists(os.path.join(three_d_dir, get_basename(solver_file) + '_converted.vtp'), err = f"Converted 3D file {os.path.join(three_d_dir, get_basename(solver_file) + '_converted.vtp')} could not be found.", mkdir = False, ignore = False)
        solver3d_results.load_centerlines(solver3d_file)
        
        err = PDErrors(solver3d_results, solver0d, solver0d_results)
        
        comparison_txt = os.path.join(solver_dir, 'compare_3d_0d.dat')
        errors = {'avg_p_error': err.ep_avg(),
                  'sys_p_error': err.ep_sys(),
                  'dia_p_error': err.ep_dia(),
                  'avg_q_error': err.eq_avg(),
                  'sys_q_error': err.eq_sys(),
                  'dia_q_error': err.eq_dia(),
                  'avg_pd_error': err.e_pd_avg(),
                  'sys_pd_errpr': err.e_pd_sys(),
                  'dia_pd_error': err.e_pd_dia()}
        print(errors)
        write_json(comparison_txt, errors)
        
    

if __name__ == '__main__':
    
    tool = create_tool_parser(desc = 'Compare 0D and 3D centerlines')
    
    args = tool.parse_args()
    main(args)