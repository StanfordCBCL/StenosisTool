
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    
    data = Path('../../data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR_old/training_results/run_32768/probability_histograms')
    
    dir1 = data / 'all_distributions' / '0'

    x = np.load(dir1 / 'data.npy', allow_pickle=True).item()
    fig, ax = plt.subplots(1,1)
    
    # Histogram
    point = 4 # mPAP
    hist, bins = np.histogram(x['yhat'][:, point] / 1333.22, bins = 'auto', density = True)
    ax.bar(bins[:-1], hist, width=np.diff(bins),edgecolor="black", align="edge")
    #plt.show()
    
    # Density Plot
    idx = np.argwhere(x['yhat'][:, point]/1333.22 > bins[0]   )
    idx = np.argwhere(x['yhat'][idx.flatten(),point]/1333.22 < bins[1])

    idxx = x['x'][idx.flatten(), 0] < .2
    print(x['x'][idxx.flatten(), :])
    exit()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x['x'][idx.flatten(),0], x['x'][idx.flatten(),1], x['x'][idx.flatten(),2])
    ax.set_xlabel("LPA Proximal")
    ax.set_ylabel("RPA Proximal")
    ax.set_zlabel("RPA Distal")
    plt.show()