# checks sanity for artificial stenosis
# results must have 1000 ts per cycle runs

from functions.solver_io import SolverResults
import sys
import re
import numpy as np


if __name__ == '__main__':
    
    results_file = sys.argv[1]
     
    results = SolverResults(results_file)
    
    for name in results.pressures:
        if re.search('BC[1-9]*_outlet', name):
            lc = -1 * 1000
            time = results.time[lc:]
            time = time - time[0] + results.time[0]
            
            pressure = np.trapz(results.pressures[name][lc:], time) / time[-1] 
            print(name, pressure / 1333.32)
            
    for name in results.flow:
        if re.search('BC[1-9]*_outlet', name):
            lc = -1 * 1000
            time = results.time[lc:]
            time = time - time[0] + results.time[0]
            
            flows = np.trapz(results.flow[name][lc:], time) / time[-1]
            print(name, flows)
        