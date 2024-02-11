
from svinterface.core.polydata import Centerlines
import numpy as np
import re
class Extract1D():
    """Extract results from 1D centerlines. """
    
    def __init__(self, centerlines: Centerlines):
        
        self.centerlines = centerlines
        
    def extract_valid(self):
        """#!INCOMPLETE"""
        
        # check its there
        names = self.centerlines.get_pointdata_arraynames()
        if 'valid' not in names:
            raise ValueError(f"'valid' not foudn in Pointdata arrays")
        
        # check valid covers junctions_0d
        valid = self.centerlines.get_pointdata_array("valid")
        gidx = self.centerlines.get_pointdata_array(self.centerlines.PointDataFields.NODEID)
        valid_gidx = gidx[np.where(valid == 1)[0]]
        
        # if it is valid
        #! considers ordered.
    
        for f in ['pressure', 'flow']:
            vals = []
            times = []
            for arr in names:
                if arr.startswith(f):
                    times.append(float(f.split('_')[-1]))
                    vals.append(self.centerlines.get_pointdata_array(arr)[valid_gidx])
        
            # integrate
            times = np.array(times)
            vals = np.array(vals)
            
            avg = None
            
        
                    
                    
            