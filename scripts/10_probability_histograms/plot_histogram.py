
#! Rewrite For plotting

# Compute Density Histogram or KDE of joint distribution
from sklearn.neighbors import KernelDensity
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc


from svinterface.utils.misc import d2m
from svinterface.plotting.params import set_params

def plot_dens(x, yhat_col, bounds, ratios=[1,1,1]):
    fig = plt.figure(figsize=(4,3.3))
    ax = plt.axes(projection='3d', computed_zorder=False)
    colors = ['b','r', 'g']
    markers = ['o','>','D']
    
    for i in range(3):
        idx1 = np.argwhere((yhat_col >= bounds[i][0]) & (yhat_col < bounds[i][1])).squeeze()
        idx1 = np.random.choice(idx1, size=int(len(idx1) * ratios[i] ), replace = False)
        ax.scatter3D(x[idx1,0], x[idx1,1], x[idx1,2], edgecolor = mplc.to_rgba('black', .05), facecolor=mplc.to_rgba(colors[i], 1) ,marker = markers[i],  label=f'[{bounds[i][0]:.2f}, {bounds[i][1]:.2f}] mmHg', zorder = 3 - i)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_xlabel("LPA Proximal")
    ax.set_ylabel("RPA Proximal")
    ax.set_zlabel("RPA Distal")
    fig.legend(fontsize = plt.rcParams['font.size'] - 2)
    ax.view_init(elev=45., azim=-45)
    
    return fig, ax

if __name__ == '__main__':
    
    # data
    # data = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/all_distributions/0/data.npy'
    # data = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/rpa_distal_conditional/0.6/data.npy'
    data = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768/probability_histograms/rpa_distal_conditional/1.0/data.npy'
    
    # # plotting parameter (data0)
    # point = 0
    # ratios = [1, 1, 1]
    # # percentage of density for overlapped bins
    # overlapped = [0, 0.001, 0.005, 0.01]
    # # diastolic
    # diastolic = True
    # bandwiths = [0.01,.1, .5]
    
    # # plotting parameter (data1)
    # point = 0
    # ratios = [.025, .025, .025]
    # # percentage of density for overlapped bins
    # overlapped = [0, 0.01, 0.05, 0.1]
    # # diastolic
    # diastolic = False
    # # bandwiths
    # bandwiths = [0.025,.05]
    
    # plotting parameter (data1)
    point = 0
    ratios = [.025, .025, .025]
    # percentage of density for overlapped bins
    overlapped = [0, 0.01, 0.05, 0.1]
    # diastolic
    diastolic = False
    # bandwiths
    bandwiths = [0.025,.05]
    
    # parse data
    data = np.load(data, allow_pickle=True).item()
    x_data = data['x']
    yhat = data['yhat']
    baseline_x = data['baseline'][0]
    baseline_yhat = data['baseline'][1][0]
    

    set_params(use_latex=True, small_ticks=True, linewidth=.1, plw = 0.01, llw = 1)
    # convert baseline to mmHg
    base_idx = np.arange(0, len(yhat[0]), 6)
    pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
    baseline_yhat[pressures_idx] = d2m(baseline_yhat[pressures_idx])
    
    yhat[:,pressures_idx] = d2m(yhat[:,pressures_idx])

    if diastolic:
        count = 3
        names = [ 'Diastolic Pressure [mmHg]','Mean Pressure [mmHg]','Systolic Pressure [mmHg]']
        overlap_start = 1
    else:
        count = 2
        names = ['Mean Pressure [mmHg]','Systolic Pressure [mmHg]']
        overlap_start = 0
            
        
    fig, ax = plt.subplots(1, count, figsize = (3*count + 1, 2.5) )
    ax = ax.flatten()

    #! when rewriting, add an option for viewing original histograms
    for j in range(count):
        point_idx =  point * 6 + 3 + j + (1-overlap_start)
        
        X = yhat[:, point_idx][:, np.newaxis]
        X_plot = np.linspace(X.min(), X.max(), 1000)[:, np.newaxis]
        # print(X_plot)
        kde = KernelDensity(kernel='gaussian', bandwidth=bandwiths[j]).fit(X)
        log_dens = kde.score_samples(X_plot)
        ax[j].axvline(x=baseline_yhat[point_idx], color = 'r', label = 'Baseline', zorder = 3)
        # ax[j].hist(X[:,0], bins = 'auto', density = True)
        x = X_plot[:, 0]
        y = np.exp(log_dens)
        y[0] = 0
        y[-1] = 0            
        ax[j].plot(x, y, linestyle='--', zorder = 2)
        ax[j].fill(x, y, alpha = 0.3, zorder = 1)
        ax[j].set_xlabel(names[j])
        ax[j].set_ylabel("Density")
        ax[j].set_ylim(bottom = 0)
        
        if j >= overlap_start:
            
   
            
            
            colors = ['b','r','g']
            bounds = []
            for k in range(3):
                cum_y = np.cumsum(y)
                cum_y_norm = cum_y / cum_y[-1]
                idxes =  np.argwhere((cum_y_norm < overlapped[k+1]) & (cum_y_norm >= overlapped[k]))
                lidx = idxes[0][0]
                uidx = idxes[-1][0]
                # print(uidx, lidx)
                new_x = np.concatenate(([x[lidx]],x[lidx:uidx + 1], [x[uidx]]) )
                new_y = np.concatenate(([0], y[lidx:uidx + 1], [0]))
                ax[j].fill(new_x, new_y, color = colors[k])
                bounds.append((x[lidx], x[uidx]))
            
            tmp_fig, tmp_ax = plot_dens(x_data, yhat[:, point_idx], bounds = bounds, ratios = ratios)
            tmp_fig.savefig(f"images/paper/temp{j}.pdf")
        
    
    ax[2 - (1-overlap_start)].legend(fontsize=plt.rcParams['font.size']-2,  loc='upper right')
    fig.tight_layout()


    # fig.savefig('images/paper/07_results/joint_dist_point_0_pressures.pdf')
    # fig.savefig('images/paper/07_results/cond_dist_point_0_pressures.0.6.pdf')
    fig.savefig('images/paper/07_results/cond_dist_point_0_pressures.1.0.pdf')
    fig.suptitle(f"Point {point}")
        
    plt.show()
   