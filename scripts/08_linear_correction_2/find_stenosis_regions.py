# Finds stenosis regions by comparing two centerlines and taking standard deviation of distances (May require several iterations of observing model)
#! a better way would just be to actually know the region of stenosis, but this is challenging


from svinterface.core.polydata import Centerlines
import numpy as np

    
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
    import argparse
    parser = argparse.ArgumentParser(description='Finds stenosis regions by comparing two centerlines and finding the top n points')
    parser.add_argument('-dir', default = "data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/3D_DIR/LPA_stent", help = '3D dir containing diseased and stented centerlines')
    parser.add_argument("-max_points", type=int, help = 'maximum number of points')
    parser.add_argument("-outdir", default='.', help='output directory')
    parser.add_argument("-name", default = 'test', help = 'file to output relevant regions to')
    
    args = parser.parse_args()
    
    name = args.name
    dir3d = Path(args.dir)
    threed_files = sorted(list((dir3d).glob("*.vtp")), key=lambda x: len(str(x)))
    stented_file = threed_files[0]
    stented_centerlines = Centerlines.load_centerlines(stented_file)
    dis_file = threed_files[2]
    dis_centerlines = Centerlines.load_centerlines(dis_file)
    
    d = get_distances(dis_centerlines, stented_centerlines)
    distances = np.array(d)[:,2]
    
    shapex = dis_centerlines.get_pointdata_array('GlobalNodeId').shape
    shapey = stented_centerlines.get_pointdata_array('GlobalNodeId').shape
    # clear arrays
    clear_other_arrays(dis_centerlines)
    clear_other_arrays(stented_centerlines)
    
    # save data
    # compute several std.
    for n in range(1,args.max_points + 1):
        # get top n distances
        t = np.array(d[:n])
        
        x = np.zeros(shapex)
        x[t[:,0].astype(int)] = 1
        dis_centerlines.add_pointdata(x, f'important_{n}')
        
        y = np.zeros(shapey)
        # print(lpa_dis_centerlines.polydata)
        y[t[:,1].astype(int)] = 1
        
        stented_centerlines.add_pointdata(y, f'important_{n}')
        
    dis_centerlines.write_polydata(str(Path(args.outdir) / f"{name}_stenosis_search_diseased.vtp"))
    stented_centerlines.write_polydata(str(Path(args.outdir) / f"{name}_stenosis_search_stented.vtp"))