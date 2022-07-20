import sys
from functions.centerlines import Centerlines
import numpy as np


def add_avg_max_mPAP(centerlines_file ):
    centerlines = Centerlines()
    centerlines.load_centerlines(centerlines_file)
    array_num = centerlines.centerlines.GetPointData().GetNumberOfArrays()
    array_names = [centerlines.centerlines.GetPointData().GetArrayName(i) for i in range(array_num)]
    mpap = {}
    maxpap = None
    cur_max = 0
    for name in array_names:
        if name.startswith('pressure'):
            p_data = centerlines.get_pointdata(name)
            if max(p_data) > cur_max:
                cur_max = max(p_data)
                maxpap = float(name.split('_')[-1])
            
            mpap[float(name.split('_')[-1])] = p_data
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([mpap[tidx]/1333.22 for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'mPAP')
    
    centerlines.add_pointdata(mpap[maxpap]/1333.22, 'maxPAP_' + str(maxpap))
    
    centerlines.write_centerlines(centerlines_file)


if __name__ == '__main__':
    
    centerline_file = sys.argv[1]
    add_avg_max_mPAP(centerline_file)
    
    
    