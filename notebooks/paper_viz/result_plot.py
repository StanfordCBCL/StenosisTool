
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    
    data = Path('../../data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms')
    
    dir1 = data / 'all_distributions' / '0'

    x = np.load(dir1 / 'data.npy', allow_pickle=True).item()
    print(x['yhat'].shape)
    fig, ax = plt.subplots(1,1)
    print("Plot")
    # ax.hist(x['yhat'][:, pressures], bins = 'auto', density=True)
    # plt.show()
    # fig.savefig('test.png')
    # fig.close()
    

    
    point = 4
    hist, bins = np.histogram(x['yhat'][:, point] / 1333.22, bins = 'auto', density = True)
    print(hist)
    ax.bar(bins[:-1], hist, width=np.diff(bins),edgecolor="black", align="edge")
    print(bins)
    plt.show()
    idx = np.argwhere(x['yhat'][:, point]/1333.22 > 31.793278   )
    print(idx)
    idx = np.argwhere(x['yhat'][idx.flatten(),point]/1333.22 < 32.001232)
    print(idx.shape)
    
    print(x['x'][idx.flatten()])
    ax = plt.axes(projection='3d')
    ax.scatter3D(x['x'][idx.flatten(),0], x['x'][idx.flatten(),1], x['x'][idx.flatten(),2])
    plt.show()