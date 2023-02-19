# File: 0D_model_to_3D.py
# File Created: Thursday, 26th January 2023 8:49:44 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 14th February 2023 1:22:39 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Takes a 0D LPN model and reconstructs in 3D what the 0D represents as a legacy vtk file
#! Untested

import argparse
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from vtk.util.numpy_support import numpy_to_vtk as n2v
from svinterface.core.zerod.lpn import LPN
from svinterface.core.polydata import Polydata, LegacyVTK



def get_ref_frame(vec):
  """
  Generate a reference frame given an axial vector vec
  """
  v1 = (vec)/np.linalg.norm(vec)
  v2 = np.cross(v1,np.array([0,0,1]))
  if(np.linalg.norm(v2) > 1.0e-4):
    v2 = np.cross(v1,np.array([0,1,0]))
  v2 = v2/np.linalg.norm(v2)
  v3 = np.cross(v1,v2)
  v3 /= np.linalg.norm(v3)
  return np.concatenate((v1.reshape(1,-1),v2.reshape(1,-1),v3.reshape(1,-1)),axis=0)


def get_points(target, origin,n_segments, n_theta, radius):
  """
  Create all ponts for the cylindrical geometry
  """
  # Get segment length
  segment_length = np.linalg.norm(target-origin)
  # Create reference frams
  # ref[0] is the axial direction 
  # ref[1] and ref[2] are the two orthogonal directions
  ref = get_ref_frame(target-origin)
  # Create points 
  points = np.zeros(((n_segments + 1)*n_theta,3))
  for loopA in range(n_segments + 1):
    curr_center = origin + loopA*(segment_length/n_segments)*ref[0]
    for loopB in range(n_theta):
      points[loopA*n_theta + loopB] = curr_center + ref[1]*radius*np.cos(2.0*np.pi*loopB/n_theta) + ref[2]*radius*np.sin(2.0*np.pi*loopB/n_theta)
  return points




if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description="Takes a 0D LPN model and reconstructs in 3D what the 0D represents")
    
    parser.add_argument("-c", dest = 'centerlines', help = 'centerlines file to figure out what the 0D is modeling.')
    parser.add_argument("-lpn", help = '0D LPN network to convert')
    parser.add_argument("-o", help = "output file destination")
    parser.add_argument("--ntheta", type = int, default = 100, help = "number of circumference points for each slice" )
    parser.add_argument("--nsegments", type = int, default = 1, help = "Number of triangular segments to split vessel lengthwise into.")
    
    args = parser.parse_args()
    
    # get lpn branch tree
    lpn = LPN.from_file(args.lpn)
    lpn_root = lpn.get_tree()
    
    # centerlines
    c = Polydata()
    c.load_polydata(args.centerlines)
    
    # new 3D model
    three_d = LegacyVTK()
    
    # load branchids
    branch_ids = c.get_pointdata_array("BranchId")
    paths = c.get_pointdata_array("Path")
    points = c.get_points()

    # iterate through each branch
    tracked_segments = 0
    for branch_node in lpn.tree_bfs_iterator(lpn_root, "branch"):
        bidx = branch_node.id
        
        # retrieve branch specific points
        branch_pointidx = np.where(branch_ids == bidx)[0]
        branch_paths = paths[branch_ids == bidx]
        
        # initialize axial vector idx's on centerlines (used to find 3D coords)
        origin_idx = None
        target_idx = branch_pointidx[0]
        branch_cur_len = 0
        for vessel in branch_node.vessel_info:
            
            # shift idx to next vessel segment in branch
            branch_cur_len += vessel['vessel_length']
            origin_idx = target_idx
            target_idx = branch_pointidx[np.where(abs(branch_paths - branch_cur_len) < 1e-8)[0]][0] # account for slight numerical inaccuracy from float comparisons

            # retrieve 3D coords of axial vector
            origin = points[origin_idx]
            target = points[target_idx]
            
            # get radius
            radius = lpn.get_vessel_radius(vessel['vessel_id'])
            
            
            # get points
            new_points = get_points(target, origin, args.nsegments, args.ntheta, radius)
            
            # get connectivity for each segment of this vessel
            connectivity = []
            n_segments = tracked_segments + args.nsegments 
            for loopA in range(tracked_segments, n_segments):
                tmp = []
                for loopB in range(args.ntheta):
                    tmp += [loopA*args.ntheta+loopB,(loopA+1)*args.ntheta+loopB]
                tmp += [loopA*args.ntheta,(loopA+1)*args.ntheta]
                connectivity.append(tmp)
            connectivity = np.array(connectivity)
            
            # add it to the vtk aggregator
            three_d.add_polydata(new_points, connectivity)
            
            tracked_segments = n_segments + 1
            
    
    three_d.write_vtk(args.o)
        
            
            
            
            
            
