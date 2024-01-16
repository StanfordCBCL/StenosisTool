import vtk
import numpy as np
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.lpn import LPN
from collections import defaultdict
from scipy.interpolate import interp1d

if __name__ == '__main__':
    base = "data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/"
    
    threed = Centerlines()
    threed.read_polydata(base + "3D_DIR/prestent/AS1_SU0308_3D_centerlines.formatted.vtp")
    
    zerod_nocor = Centerlines()
    zerod_nocor.read_polydata(base + "LPN_DIR/AS1_SU0308.sim.0/centerline_projection.vtp")
    
    zerod_final = Centerlines()
    zerod_final.read_polydata(base + "LPN_DIR/AS1_SU0308.sim.5/centerline_projection.vtp")
    
    err = Centerlines()
    err.read_polydata(base + "3D_DIR/prestent/AS1_SU0308_3D_prestent_centerlines.full.formatted.vtp")
    for name in err.get_pointdata_arraynames():
        err.remove_pointdata_array(name)
    
    # lpn
    lpn = LPN.from_file("data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/LPN_DIR/AS1_SU0308.sim.5/AS1_SU0308.in")
    
    # err.write_polydata("images/plot_data/err_plot.vtp")
    
    ##! Get results in python format
    branch_results = {}    
    # write time
    branch_results['time'] = ['err1', 'err2']
    
    branch_results['distance'] = {}
    branch_results['errors'] = {}
    
    branch_tree = lpn.get_tree()
    
    threed_ap = threed.get_pointdata_array('avg_pressure')
    zerod_final_ap = zerod_final.get_pointdata_array('avg_pressure')
    zerod_nocor_ap = zerod_nocor.get_pointdata_array('avg_pressure')
    
    for node in lpn.tree_bfs_iterator(branch_tree, allow='branch'):
        branch_id = node.id
        
        # gids    
        gids = node.vessel_info[0]['gid']

        # distances
        vess_D = [0.0]
        vess_D.append(node.vessel_info[0]['vessel_length'])
        
        if len(node.ids) > 1:
            for idx, v in enumerate(node.vessel_info[1:], 1):
                gids.append(v['gid'][1])
                vess_D.append(vess_D[-1] + v['vessel_length'])
        
        # get initial values for first vessel
        vess_E = [[np.abs(threed_ap[g] - zerod_nocor_ap[g]), np.abs(threed_ap[g] - zerod_final_ap[g])] for g in gids]

    
        branch_results['distance'][branch_id] = np.array(vess_D)
        branch_results['errors'][branch_id] = np.array(vess_E)
   
    results = branch_results
    
    
    
    # # assemble output dict
    rec_dd = lambda: defaultdict(rec_dd)
    arrays = rec_dd()

    # # extract point arrays from geometries
    arrays_cent = {}
    for arr_name in threed.get_pointdata_arraynames():
        arrays_cent[arr_name] = threed.get_pointdata_array(arr_name)

    # # centerline points
    points = err.get_points()

    # # all branch ids in centerline
    ids_cent = np.unique(arrays_cent['BranchId']).tolist()
    ids_cent.remove(-1)
    
    # # retrieve time
    times = results['time']

    # # loop all result fields
    for f in ['errors']:
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
        
    # add arrays to centerline and write to file
    for f, a in arrays.items():
        err.add_pointdata(array=a, array_name=f)
    
    err.write_polydata("images/plot_data/err_plot.vtp")