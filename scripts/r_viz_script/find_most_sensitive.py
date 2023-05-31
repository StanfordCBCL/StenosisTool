
import sys
import numpy as np

from svinterface.manager.baseManager import Manager
from svinterface.core.zerod.lpn import LPN
from svinterface.core.zerod.solver import Solver0Dcpp


if __name__ == '__main__':
    
    M = Manager(sys.argv[1])
    
    lpn = LPN.from_file(M['workspace']['lpn'])
    
    s = Solver0Dcpp(lpn, mean_only=True)
    
    # initial
    init_results = s.run_sim()

    # get  tree to loop through each
    local_flow_change = np.zeros(lpn.num_vessels())
    local_pressure_change = np.zeros(lpn.num_vessels())
    tree = lpn.get_tree()
    for vess in lpn.tree_bfs_iterator(tree, allow = 'branch'):
        print("Running", vess.ids)
        for idx, vidx in enumerate(vess.ids):
            # half the R
            print(vess.vessel_info[idx]["zero_d_element_values"]['R_poiseuille'], vess.vessel_info[idx]["zero_d_element_values"]['R_poiseuille'] * .5, )
            lpn.change_vessel(vessel_id_or_name=vidx, R = vess.vessel_info[idx]["zero_d_element_values"]['R_poiseuille'] * 2)
            
        s = Solver0Dcpp(lpn, mean_only=True)
        results = s.run_sim()
        
        # reverse half by doubling
        for idx, vidx in enumerate(vess.ids):
            local_pressure_change[vidx] = results.result_df.iloc[vidx]['pressure_out'] - init_results.result_df.iloc[vidx]['pressure_out']
            local_flow_change[vidx] = results.result_df.iloc[vidx]['flow_out'] - init_results.result_df.iloc[vidx]['flow_out']
            lpn.change_vessel(vessel_id_or_name=vidx, R = vess.vessel_info[idx]["zero_d_element_values"]['R_poiseuille'] * .5)
            
    
    np.save('flow_test2.npy', local_flow_change)
    np.save('pressure_test2.npy', local_pressure_change)
    init_results.save_csv('init_res2.csv')