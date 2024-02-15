# Compute Density Histogram or KDE of joint distribution
from sklearn.neighbors import KernelDensity
import numpy as np
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

from sample_data import RepairDistribution

from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Plot Density Histogram")
    parser.add_argument("-data", default='data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy', help='path to data.npy generated from sampling')
    parser.add_argument("-outdir", default='.', help='directory to output results to (will create if it doesnt exist)')
    parser.add_argument("-prefix", default='', help='prefix to add to file output (used to label)')
    parser.add_argument("-kde", action='store_true', help='use the kde rather than histogram (kde will probably have issues with bandwith)')
    parser.add_argument("-points", nargs='+', default=[0], type=int, help='number of plots to limit to, rather than plotting entire set of points')
    parser.add_argument("--p", action='store_true', help='Plot pressures')
    parser.add_argument("--q", action='store_true', help='Plot flows')
    parser.add_argument("--latex", action='store_true', help='Plot using Latex (nicer but needs Latex installed)')
    
    args = parser.parse_args()
    # parse data
    data = np.load(args.data, allow_pickle=True).item()
    x = data['x']
    yhat = data['yhat']
    baseline_x = data['baseline'][0]
    baseline_yhat = data['baseline'][1][0]
    
    # suffic for saving
    if args.p and args.q:
        suffix = ''
    elif args.p:
        suffix = '_pressures'
    elif args.q:
        suffix = '_flows'   
        
    # prefix
    prefix = args.prefix if args.prefix == '' else args.prefix + '_'


    # Density Hist
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    set_params(use_latex=args.latex, small_ticks=True, linewidth=.1)
    if not args.kde:
        # convert baseline to mmHg
        base_idx = np.arange(0, len(yhat[0]), 6)
        pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
        baseline_yhat[pressures_idx] = d2m(baseline_yhat[pressures_idx])
        
        # density hist    
        dist = RepairDistribution(3)
        hist = dist.get_histograms(yhat, points=args.points)
        
        for pidx, point in enumerate(args.points):
            fig = dist.plot_single_histogram(hist[pidx], p=args.p, q=args.q, baseline_yhat=baseline_yhat[point*6:point*6+6])  
            fig.savefig(str(outdir / f'{prefix}point_{point}{suffix}.png'))
            fig.suptitle(f"Point {point}")
            
    else:
        # kde
        #convert pressure data to mmHg
        base_idx = np.arange(0, len(yhat[0]), 6)
        pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
        yhat[:, pressures_idx] = d2m(yhat[:, pressures_idx])
        
        p = args.p; q = args.q
        for point in args.points:
            if p and q:
                fig, ax = plt.subplots(2, 3, figsize = (10, 5) )
                ax = ax.flatten()
                names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)', 'Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
                bandwiths = [.01, .1, .3, .01, .1, .3]
                for j in range(6):
                    X = yhat[:, point + j][:, np.newaxis]
                    X_plot = np.linspace(X.min(), X.max(), 1000)[:, np.newaxis]
                    kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[j]).fit(X)
                    log_dens = kde.score_samples(X_plot)
                    x, y = X_plot[:, 0], np.exp(log_dens)
                    y[0] = 0
                    y[-1] = 0
                    ax[j].plot(x,y)
                    ax[j].set_xlabel(names[j])
                    ax[j].set_ylabel("Density")
                fig.tight_layout()
                
            elif p:
                fig, ax = plt.subplots(1, 3, figsize = (10, 2.5) )
                ax = ax.flatten()
                names = ['Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
                bandwiths = [.01, .1, .3]
                for j in range(3, 3+len(names)):
                    X = yhat[:, point + j][:, np.newaxis]
                    X_plot = np.linspace(X.min(), X.max(), 1000)[:, np.newaxis]
                    kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[j-3]).fit(X)
                    log_dens = kde.score_samples(X_plot)
                    x, y = X_plot[:, 0], np.exp(log_dens)
                    y[0] = 0
                    y[-1] = 0
                    ax[j-3].plot(x,y)
                    ax[j-3].set_xlabel(names[j-3])
                    ax[j-3].set_ylabel("Density")
                fig.tight_layout()

            elif q:
                fig, ax = plt.subplots(1, 3, figsize = (10, 2.5) )
                ax = ax.flatten()
                names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)']
                bandwiths = [.01, .1, .3]
                for j in range(len(names)):
                    X = yhat[:, point + j][:, np.newaxis]
                    X_plot = np.linspace(X.min(), X.max(), 1000)[:, np.newaxis]
                    kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[j]).fit(X)
                    log_dens = kde.score_samples(X_plot)
                    x, y = X_plot[:, 0], np.exp(log_dens)
                    y[0] = 0
                    y[-1] = 0
                    ax[j].plot(x,y)
                    ax[j].set_xlabel(names[j])
                    ax[j].set_ylabel("Density")
                fig.tight_layout()

    plt.show()
   