# File: centerline_match.py
# File Created: Wednesday, 14th June 2023 5:11:09 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:24:41 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Determines corresponding GID on the stented model centerlines that are relevant to the prestented model centerlines


from svinterface.manager import Manager
from svinterface.core.polydata import Centerlines

import argparse
import numpy as np
from tqdm import tqdm

def dist(p1, p2):
    # euclidian distance
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2 )
    
def find_closest(p1, stented_gids, stented_points):
    """Finds GID of closest point"""
    min_gid = stented_gids[0]
    min_dist = dist(p1, stented_points[min_gid])
    for gid, p2 in list(zip(stented_gids, stented_points)):
        tmp_dist = dist(p1, p2)
        # check if smaller. update gid
        if tmp_dist < min_dist:
            min_dist = tmp_dist
            min_gid = gid
    
    return min_gid
            
        

def match_centerlines(diseased_cent: Centerlines, stented_cent: Centerlines):
    # check for valid
    caps = diseased_cent.get_pointdata_array("Caps_0D")
    junc = diseased_cent.get_pointdata_array("Junctions_0D")
    vess = diseased_cent.get_pointdata_array("Vessels_0D")
    valid = caps + junc + vess + 3
    
    # get valid gid
    gids = diseased_cent.get_pointdata_array("GlobalNodeId")
    valid_gids = gids[valid > 0]
    poi_points = diseased_cent.get_points()[valid_gids]
    
    # for every point on the stented cent, find closest point
    stented_points = stented_cent.get_points()
    stented_gids = stented_cent.get_pointdata_array("GlobalNodeId")
    
    gids_map = {} # old_gid: new_gid
    for idx, p1 in tqdm(enumerate(poi_points), total = len(poi_points)):
        gids_map[valid_gids[idx]] = find_closest(p1, stented_gids, stented_points)
    
    # map gids func
    gid_mapper = lambda new_gid: gids_map[new_gid]
    
    # decompose into caps junc and vess
    mapped_caps = np.zeros_like(stented_gids) - 1
    caps_gids = gids[caps > -1]
    mapped_caps[np.array([gid_mapper(g) for g in caps_gids])] = caps_gids
    stented_cent.add_pointdata(mapped_caps, 'Caps_0D')
    
    mapped_juncs = np.zeros_like(stented_gids) - 1
    juncs_gids = gids[junc > -1]
    mapped_juncs[np.array([gid_mapper(g) for g in juncs_gids])] = juncs_gids
    stented_cent.add_pointdata(mapped_juncs, 'Junctions_0D')
    
    mapped_vess = np.zeros_like(stented_gids) - 1
    vess_gids = gids[vess > -1]
    mapped_vess[np.array([gid_mapper(g) for g in vess_gids])] = vess_gids
    stented_cent.add_pointdata(mapped_vess, 'Vessels_0D')
    
    return
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Determines the stented model centerline GID that corresponds to relevant GID in diseased model")
    parser.add_argument("-i", help = "Config file")
    parser.add_argument("-c", help = "path to stented centerlines")
    args = parser.parse_args()
    
    M = Manager(args.i)
    
    #load centerlines
    dis_cent = Centerlines.load_centerlines(M['workspace']['centerlines'])
    stent_cent = Centerlines.load_centerlines(args.c)
    
    # match
    match_centerlines(diseased_cent=dis_cent,
                      stented_cent=stent_cent)
    
    # update
    stent_cent.write_polydata(args.c)
    
    
    
    
    
    