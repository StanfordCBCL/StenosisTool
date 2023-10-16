# Plot 3D density plots (only useful for N=3. For larger values, this will be less useful)

import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mp
from sklearn.neighbors import KernelDensity
from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params

if __name__ == '__main__':
    point = 3
    dfiles = ['data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy',
              'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/1/data.npy',
              'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/2/data.npy']
    set_params(linewidth=.1, use_latex=True, small_ticks=True, )
    fig, ax = plt.subplots(2, 3, figsize = (10, 5) )
    ax = ax.flatten()
    first = True
    counter = -1
    labels = [ '$\sigma=0$', '$\sigma=1$', '$\sigma=2$']
    marker = ['--', '-.', ':']
    bandwiths = [ [.01, .5, 1, .01, .4, .6], [.1, .4, 1, .1, .4, .6], [.2, .6, 1, .2, .4, .6]]
    for datafile in dfiles:
        counter += 1
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
        
        
        

        names = ['Diastolic RPA Flow (mL/s)', 'Mean RPA Flow (mL/s)', 'Systolic RPA Flow (mL/s)', 'Diastolic PAP (mmHg)','Mean PAP (mmHg)','Systolic PAP (mmHg)']
        for j in range(6):
            if j < 3:
                point = 82
                h = hist[point]
                
            else:
                point = 0
                h = hist[point]
                
            X = yhat[:, point * 6 + j][:, np.newaxis]
            X_plot = np.linspace(X.min(), X.max(), 100)[:, np.newaxis]
            kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[counter][j]).fit(X)
            log_dens = kde.score_samples(X_plot)
            
            if first:
                ax[j].axvline(x=baseline_yhat[point * 6 + j], color = 'r', label = 'Baseline', zorder = 3)
            
            x = X_plot[:, 0]
            y = np.exp(log_dens)
            y[0] = 0
            y[-1] = 0            
            ax[j].plot(x, y, label = labels[counter], linestyle=marker[counter], zorder = 2)
            ax[j].fill(x, y, alpha = 0.3, zorder = 1)
        
            
            #ax[j].bar(h[j][1][:-1], h[j][0], width=np.diff(h[j][1]), alpha = 1, align="edge")

            ax[j].set_xlabel(names[j])
            ax[j].set_ylabel("Density")
            ax[j].set_ylim(bottom = 0)
        first = False
    ax[0].set_ylim(0, .5)
    ax[3].set_ylim(0, .5)
    ax[2].legend(fontsize=plt.rcParams['font.size']-2)
    mp.rcParams['axes.linewidth'] = .1
    mp.rcParams['lines.linewidth'] = .1
    mp.rcParams['patch.linewidth'] = .1
    fig.tight_layout()
    
    fig.savefig("images/paper/07_results/joint_dist_all_std.pdf")
        
    plt.show()