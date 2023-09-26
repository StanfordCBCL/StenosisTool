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