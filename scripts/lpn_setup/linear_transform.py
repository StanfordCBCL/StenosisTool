# File: linear_transform.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Saturday, 18th February 2023 5:11:05 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear transform on the junctions






import argparse

from svinterface.core.zerod.lpn import LPN
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
import numpy as np
from concurrent.futures import ProcessPoolExecutor


def linear_transform(zerod_lpn: LPN, threed_c: Centerlines):
    
    # solve for original solver
    original = Solver0Dcpp(zerod_lpn, last_cycle_only=True,mean_only=True, debug = True)
    original_sim = original.run_sim()
    
    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect
    junction_outlet_vessels = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_outlet_vessels += junc_node.vessel_info[0]['outlet_vessels']

    # get pressures
    pressures = original_sim.result_df.iloc[junction_outlet_vessels]['pressure_in'].as_numpy()

    
    
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Perform a linear optimization on the branches")
    parser.add_argument("-3d", dest = 'threed', help = '3d centerlines file')
    parser.add_argument("-0d", dest = 'zerod', help = "0d LPN")
    
    args = parser.parse_args()
    
    
    # get LPN and convert to BVJ
    zerod_lpn = LPN.from_file(args.zerod)
    zerod_lpn.normal_to_bvj()
    
    # load centerlines
    threed_c = Centerlines.load_centerlines(args.threed)
    
    linear_transform(zerod_lpn,threed_c)