from svinterface.core.polydata import Centerlines

import numpy as np


def map_stented_to_prestent(stented: Centerlines, prestent: Centerlines):
    """ Performs map of stented to prestented"""
    j = stented.get_pointdata_array('Junctions_0D')
    v = stented.get_pointdata_array('Vessels_0D')
    c = stented.get_pointdata_array('Caps_0D')
    
    def add_valid(p_arr, s_arr):
        """Find valid regions"""
        jidx = np.where(j > -1)
        p_arr[j[jidx]] = s_arr[jidx]
        vidx = np.where(v > -1)
        p_arr[v[vidx]] = s_arr[vidx]
        cidx = np.where(c > -1)
        p_arr[c[cidx]] = s_arr[cidx]
        return p_arr
    
    
def get_distances(diseased_cent: Centerlines, stented_cent: Centerlines):
    '''Retrieves the distance between diseased and stented centerlines'''
    # check for valid
    caps = stented_cent.get_pointdata_array("Caps_0D")
    junc = stented_cent.get_pointdata_array("Junctions_0D")
    vess = stented_cent.get_pointdata_array("Vessels_0D")
    caps_idx = np.where(caps>-1)[0]
    junc_idx = np.where(junc >-1)[0]
    vess_idx = np.where(vess >-1)[0]
    # get all valid indices, ensuring ordering is maintained correctly even after unique sorting
    stent_indices = np.unique(np.concatenate([caps_idx, junc_idx, vess_idx]), return_inverse=True)
    reconstruct_stent = stent_indices[1]
    stent_indices = stent_indices[0]
    dis_indices = np.zeros_like(stent_indices)
    dis_indices[reconstruct_stent] = np.concatenate([caps[caps_idx],junc[junc_idx],vess[vess_idx]])
    assert len(stent_indices) == len(dis_indices), 'Length of stented and disease indices do not match up.'
    
    # get poi on both
    stented_poi = stented_cent.get_points()[stent_indices]
    dis_poi = diseased_cent.get_points()[dis_indices]
    
    def dist(p1, p2):
        # euclidian distance
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2 )

    # for every point on the stented cent, find closest point
    distances = []
    for idx, (p1, p2) in enumerate(zip(stented_poi, dis_poi)):
        distances.append((dis_indices[idx], stent_indices[idx], dist(p1, p2)))
        
    return sorted(distances, key=lambda x: x[2], reverse=True)

def clear_other_arrays(c: Centerlines):
    for name in c.get_pointdata_arraynames():
        c.remove_pointdata_array(name)
        
if __name__ == '__main__':
    from pathlib import Path
    
    dir3d = Path("data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/3D_DIR/")
    for name, std in list(zip(['LPA_stent', 'RPA_stent', 'RPA_2_stent'], [
        1, 2.5, 2])):
        threed_files = sorted(list((dir3d/name).glob("*.vtp")), key=lambda x: len(str(x)))
        lpa_stented_file = threed_files[0]
        lpa_stented_centerlines = Centerlines.load_centerlines(lpa_stented_file)
        lpa_file = threed_files[1]
        lpa_dis_centerlines = Centerlines.load_centerlines(lpa_file)
        
        d = get_distances(lpa_dis_centerlines, lpa_stented_centerlines)
        distances = np.array(d)[:,2]
        t = np.array(d[:np.where(distances > distances.mean() + std*distances.std())[0][-1] + 1])

        
        x = np.zeros_like(lpa_dis_centerlines.get_pointdata_array('GlobalNodeId'))
        x[t[:,0].astype(int)] = 1
        clear_other_arrays(lpa_dis_centerlines)
        lpa_dis_centerlines.add_pointdata(x, 'important')
        lpa_dis_centerlines.write_polydata(f"images/plot_data/{name}_diseased.vtp")
        
        x = np.zeros_like(lpa_stented_centerlines.get_pointdata_array('GlobalNodeId'))
        # print(lpa_dis_centerlines.polydata)
        x[t[:,1].astype(int)] = 1
        clear_other_arrays(lpa_stented_centerlines)
        lpa_stented_centerlines.add_pointdata(x, 'important')
        lpa_stented_centerlines.write_polydata(f"images/plot_data/{name}_stented.vtp")