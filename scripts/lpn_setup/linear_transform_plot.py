import sys
import numpy as np
import matplotlib.pyplot as plt



if __name__ == '__main__':
    
    tfile = sys.argv[1]
    data = np.load(tfile, allow_pickle=True)
    
    for key, value in data.items():
        print(key, value)