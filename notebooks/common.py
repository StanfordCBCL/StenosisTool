import numpy as np
from collections import defaultdict

def get_dict(fpath):
    """
    Read .npy dictionary saved with numpy.save if path is defined and exists
    Args:
        fpath: path to .npy file
    Returns:
        dictionary
    """
    if fpath is not None and os.path.exists(fpath):
        return np.load(fpath, allow_pickle=True).item()
    else:
        return {}
    

def rec_dict():
    """
    Recursive defaultdict
    """
    return defaultdict(rec_dict)