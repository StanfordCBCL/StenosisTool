from scipy.misc import central_diff_weights
from src.solver import Solver0D
from src.solver_results import SolverResults
from src.centerlines import Centerlines
from src.misc import get_basename, d2m
import numpy as np
import argparse
import os


def get_all_segments(solver:Solver0D, centerlines: Centerlines):
    # TODO: Use the RCR map in the solver file, iterate through the tree, and obtain the branch numbers of the caps as well as their vessel ID. Then take the vtp file for the 3D and identify the last point in the branch, Return a list of their two locations, in pairs
    outlet_branch_ids = []
    outlet_vessels = []
    tree = solver.get_branch_tree()
    for node in solver.tree_bfs_iterator(tree):
        if 'boundary_conditions' in node.vessel_info[-1]:
            outlet_branch_ids.append(node.branch_id)
            outlet_vessels.append(node.vess_id[-1])
            
    outlet_idx = []
    branchid = centerlines.get_pointdata(centerlines.PointDataFields.BRANCHID)
    for br in outlet_branch_ids:
        outlet_idx.append(np.where(branchid == br)[0][-1])
    inlet_idx = 0
    inlet_vessel = 0
    return inlet_idx, outlet_idx, inlet_vessel, outlet_vessels

def get_summary(vessdf, col):
    time = np.array(vessdf['time'])
    meanp = np.trapz(time, np.array(vessdf[col])) / (time[-1] - time[0])
    maxp = vessdf[col].max()
    minp = vessdf[col].min()
    return maxp, meanp,minp

def compute_errors(solver: Solver0D, solver_results: SolverResults, centerlines: Centerlines):
    #TODO: Retrieving pressure/flow values at particular points in particular arrays to compare values at caps
    
    inlet_idx, outlet_idx, inlet_vessel, outlet_vessels = get_all_segments(solver, centerlines)
    # find array names for max and min in 3D
    array_names = centerlines.get_pointdata_arraynames()
    for n in array_names:
        if 'maxPAP' in n:
            max_array_name = n
        if 'minPAP' in n:
            min_array_name = n
            

    all_caps = np.array([inlet_idx] + outlet_idx )
            
    max3D = centerlines.get_pointdata(max_array_name)[all_caps]
    mean3D = centerlines.get_pointdata('mPAP')[all_caps]
    min3D = centerlines.get_pointdata(min_array_name)[all_caps]
    
    

    # compute inlet
    v0 = solver_results.only_last_cycle(solver_results.vessel_df('V0'), solver.inflow.tc)
    max0, mean0, min0 =  get_summary(v0, 'pressure_in')
    
    max0D = [max0]
    mean0D = [mean0]
    min0D = [min0]
    for vidx in outlet_vessels:
        vt = solver_results.only_last_cycle(solver_results.vessel_df('V' + str(vidx)), solver.inflow.tc)
        maxt, meant, mint =  get_summary(vt, 'pressure_out')
        max0D.append(maxt)
        mean0D.append(meant)
        min0D.append(mint)
    
    max0D = d2m(np.array(max0D))
    mean0D = d2m(np.array(mean0D))
    min0D = d2m(np.array(min0D))
    
    assert len(max3D) == len(max0D)
    assert len(mean3D) == len(mean0D)
    assert len(min3D) == len(min0D)

    
    err_max = (np.abs(max0D - max3D).sum() / max3D.sum()) / len(max3D)
    err_min = (np.abs(min0D - min3D).sum() / min3D.sum()) / len(max3D)
    err_avg = (np.abs(mean0D - mean3D).sum() / mean3D.sum()) / len(max3D)
    
    
    
    
    max_0D_dp = max0D[0] - max0D[1:]
    max_3D_dp = max3D[0] - max3D[1:]
    err_max_dp = (np.abs(max_0D_dp - max_3D_dp).sum() / max_3D_dp.sum()) / len(max0D[1:])
    mean_0D_dp = mean0D[0] - mean0D[1:]
    mean_3D_dp = mean3D[0] - mean3D[1:]
    err_mean_dp = (np.abs(mean_0D_dp - mean_3D_dp).sum() / mean_3D_dp.sum()) /len(mean0D[1:])
    min_0D_dp = min0D[0] - min0D[1:]
    min_3D_dp = min3D[0] - min3D[1:]
    err_min_dp = (np.abs(min_0D_dp - min_3D_dp).sum() / min_3D_dp.sum())/len(min0D[1:])
    
    
    

    
    
    
    
    
    return err_max, err_avg, err_min, err_max_dp, err_mean_dp, err_min_dp
    
    
    
    
    

def main(args):
    # files
    solver_file = args.solver_file
    solver_result_csv = os.path.join(os.path.dirname(solver_file), get_basename(solver_file) + '_branch_results.csv')
    three_d = args.three_d
    
    solver0D = Solver0D()
    solver0D.read_solver_file(solver_file)
    
    results = SolverResults.load_from_csv(solver_result_csv)
    
    c = Centerlines()
    c.load_centerlines(three_d)
    
    err_max, err_avg, err_min, err_max_dp, err_mean_dp, err_min_dp = compute_errors(solver0D,results, c)
    print(f'Err_Sys: {err_max} Err_Dia: {err_min} Err_avg: {err_avg}')
    print(f'Err_Sys_PD: {err_max_dp} Err_Dia_PD: {err_min_dp} Err_avg_PD: {err_mean_dp}')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest = 'solver_file', help = '0D solver file containing vasculature to compute error for')
    parser.add_argument('-t', dest = 'three_d', help= '1D mapped representation of the 3D model')

    args = parser.parse_args()
    
    main(args)