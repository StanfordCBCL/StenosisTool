# Plot 3D density plots (only useful for N=3. For larger values, this will be less useful)

import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mp
from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params

if __name__ == '__main__':
    point = 3
    dfiles = ['data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy',
              'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/1/data.npy',
              'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/3/data.npy']
    set_params(linewidth=.1, use_latex=True, small_ticks=True, )
    fig, ax = plt.subplots(2, 3, figsize = (10, 5) )
    ax = ax.flatten()
    for datafile in dfiles:
        
        data = np.load(datafile, allow_pickle=True).item()
        x = data['x']
        yhat = data['yhat']
        #convert pressure data to mmHg
        base_idx = np.arange(0, len(yhat[0]), 6)
        pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
        yhat[:, pressures_idx] = d2m(yhat[:, pressures_idx])
        
        baseline_x = data['baseline'][0]
        baseline_yhat = data['baseline'][1][0]
        baseline_yhat[pressures_idx] = d2m(baseline_yhat[pressures_idx])
        
        hist = [[np.histogram(yhat[:, i + j], bins = 'auto', density=True) for j in range(3)]+[np.histogram(yhat[:, i + j], bins = 'auto', density=True) for j in range(3,6)] for i in range(0, len(yhat[0]), 6)]
        
        h = hist[point]
        
        

        names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)', 'Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
        for j in range(6):
            ax[j].bar(h[j][1][:-1], h[j][0], width=np.diff(h[j][1]),edgecolor="black", linewidth=.1, align="edge")
            ax[j].set_xlabel(names[j])
            ax[j].set_ylabel("Density")
    ax[0].set_ylim(0, .5)
    ax[3].set_ylim(0, .5)
    fig.legend(labels=['0 STD', '1 STD', '3 STD'], fontsize=plt.rcParams['font.size']-2)
    mp.rcParams['axes.linewidth'] = .1
    mp.rcParams['lines.linewidth'] = .1
    mp.rcParams['patch.linewidth'] = .1
    fig.tight_layout()
    
    fig.savefig("images/paper/07_results/joint_dist_all_std.png")
        
    plt.show()