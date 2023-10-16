# Plot 3D density plots (only useful for N=3. For larger values, this will be less useful)

import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params
from sklearn.neighbors import KernelDensity

def plot_dens(yhat, point, hist, meas, bins, rat):
    fig = plt.figure(figsize=(4,3.3))
    ax = plt.axes(projection='3d')
    point_idx = point*6 + meas
    # mPAP
    yhat_col = np.array(yhat[:,])
    X = yhat[:, point_idx][:, np.newaxis]
    X_plot = np.linspace(X.min(), X.max(), 100)[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[j]).fit(X)
    log_dens = kde.score_samples(X_plot)
    ax[j].axvline(x=baseline_yhat[point_idx], color = 'r', label = 'Baseline', zorder = 3)
    x = X_plot[:, 0]
    y = np.exp(log_dens)
   
   
    colors = ['b','r', 'm']
    markers = ['o','>','D']
    
    
    
    
    for i in range(len(bins) - 1):
        idx1 = np.argwhere((yhat_col > hist[point][meas][1][bins[i]]) & (yhat_col < hist[point][meas][1][bins[i + 1]])).squeeze()
        idx1 = np.random.choice(idx1, size=int(len(idx1) * rat ))
        ax.scatter3D(x[idx1,0], x[idx1,1], x[idx1,2], alpha = .5, color=colors[i],marker = markers[i],edgecolor='black', linewidth = 0.01, label=f'[{hist[point][meas][1][bins[i]]:.2f}, {hist[point][meas][1][bins[i+1]]:.2f}] mmHg')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_xlabel("LPA Proximal")
    ax.set_ylabel("RPA Proximal")
    ax.set_zlabel("RPA Distal")
    fig.legend(fontsize = plt.rcParams['font.size'] - 2)
    ax.view_init(elev=45., azim=-45)
    return fig, ax
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Plot density")
    parser.add_argument("-data", default='data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy', help='data.npy file to plot')
    parser.add_argument("-ofile", default=None, help='Outfile')
    parser.add_argument("-perc", type=float, nargs = 2, default = [0, .1], help='the density region to use for analysis')
    parser.add_argument("-point", default=0, type=int, help='Point to take density plot at')
    parser.add_argument("-rat_points", default=.025, type=float, help='ratio to reduce the number of plotted points by, if the plot is too dense. Randomly Samples.')
    meas = parser.add_mutually_exclusive_group()
    meas.add_argument("-dp", action='store_true', default=False, help='Diastolic Pressures')
    meas.add_argument("-mp", action='store_true', default=False, help='Mean Pressures')
    meas.add_argument("-sp", action='store_true', default=False, help='Systolic Pressures')
    meas.add_argument("-df", action='store_true', default=False, help='Diastolic Flows')
    meas.add_argument("-mf", action='store_true', default=False, help='Mean Flows')
    meas.add_argument("-sf", action='store_true', default=False, help='Systolic Flows')
    
    args = parser.parse_args()
    
    datafile = args.data
    
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
    point = args.point # MPA
    if args.df:
        meas = 0
    elif args.mf:
        meas = 1
    elif args.sf:
        meas = 2
    elif args.dp:
        meas = 3
    elif args.mp:
        meas = 4
    elif args.sp:
        meas = 5
    else:
        raise Exception("Invalid Measurement. (pick d, m, s flows or pressures)")
    fig, ax = plot_dens(yhat, point, hist, meas, args.perc, args.rat_points)
    ax.view_init(elev=45., azim=-45)
    if args.ofile:
        fig.savefig(args.ofile)
    plt.set_
    
    plt.show()