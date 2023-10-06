# Compute Density Histogram or KDE of joint distribution
# Copy over to 10_probability_histograms to plot
from sklearn.neighbors import KernelDensity
import numpy as np
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

from sample_data import RepairDistribution

from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params


if __name__ == '__main__':
    
    # data
    data = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy'
    # data = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/rpa_distal_conditional/0.6/data.npy'
    
    # plotting parameter
    point = 0
    p = True
    q = False
    # for overlapped bins
    interval = 1
    
    # parse data
    data = np.load(data, allow_pickle=True).item()
    x = data['x']
    yhat = data['yhat']
    baseline_x = data['baseline'][0]
    baseline_yhat = data['baseline'][1][0]
    

    set_params(use_latex=True, small_ticks=True, linewidth=.1)
    # convert baseline to mmHg
    base_idx = np.arange(0, len(yhat[0]), 6)
    pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
    baseline_yhat[pressures_idx] = d2m(baseline_yhat[pressures_idx])

    # density hist    
    dist = RepairDistribution(3)
    hist = dist.get_histograms(yhat, points=[point])

    fig = dist.plot_single_histogram(hist[0], p=p, q=q, baseline_yhat=baseline_yhat[point*6:point*6+6])  
    ax = fig.get_axes()
    colors = ['b','r', 'm']

    for j in range(1,3):
        for k in range(3):
            ax[j].bar(hist[0][j + 3][1][interval * k: interval * (k + 1) ],  hist[0][j+3][0][interval * k: interval * (k + 1)], width=np.diff(hist[0][j+3][1][interval * k: interval * (k + 1) + 1]),color = colors[k],edgecolor="black", linewidth=.1, align="edge")
    
    fig.savefig('images/paper/07_results/joint_dist_point_0_pressures.png')
    fig.suptitle(f"Point {point}")
        
    plt.show()
   