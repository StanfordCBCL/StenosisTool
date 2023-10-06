# Plot 3D density plots (only useful for N=3. For larger values, this will be less useful)

import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params

def plot_dens(yhat, point, meas, bins):
    fig = plt.figure(figsize=(4,3.3))
    ax = plt.axes(projection='3d')
    # mPAP
    yhat_col = yhat[:, point*6 + meas]
    for i in bins:
        idx1 = np.argwhere((yhat_col > hist[point][meas][1][i]) & (yhat_col < hist[point][meas][1][i+1]))[0]
        ax.scatter3D(x[idx1,0], x[idx1,1], x[idx1,2], alpha = .5, label=f'[{hist[point][meas][1][i]:.2f}, {hist[point][meas][1][i+1]:.2f}] mmHg')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_xlabel("LPA Proximal")
    ax.set_ylabel("RPA Proximal")
    ax.set_zlabel("RPA Distal")
    fig.legend(fontsize = plt.rcParams['font.size'] - 2)

    return fig, ax
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Plot density")
    parser.add_argument("-data", default='data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy', help='data.npy file to plot')
    parser.add_argument("-range")
    datafile = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy'
    
    data = np.load(datafile, allow_pickle=True).item()
    x = data['x']
    yhat = data['yhat']
    #convert pressure data to mmHg
    base_idx = np.arange(0, len(yhat[0]), 6)
    pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
    yhat[:, pressures_idx] = d2m(yhat[:, pressures_idx])
    
    baseline_x = data['baseline'][0]
    baseline_yhat = data['baseline'][1][0]
    
    hist = [[np.histogram(yhat[:, i + j], bins = 'auto', density=True) for j in range(3)]+[np.histogram(yhat[:, i + j], bins = 'auto', density=True) for j in range(3,6)] for i in range(0, len(yhat[0]), 6)]
    
    
    set_params(use_latex=True, small_ticks=True, )
    
    # Density Plot
    point = 0 # MPA
    meas = 5
    bins = range(20,23)
    fig, ax = plot_dens(yhat, point, meas, bins)
    ax.view_init(elev=30., azim=-45)
    #fig.savefig("images/paper/07_results/joint_dist_0std_midbins_sPAP.pdf")
    
    
    plt.show()