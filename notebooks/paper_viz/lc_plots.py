import matplotlib.pyplot as plt
import numpy as np
import matplotlib

fs=12
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('text', usetex=True)

def plot_valid2(threed_data, zerod_data):
    
    results_3d = np.load(threed_data, allow_pickle=True).item()
    threed_mins = results_3d['3d_mins']
    threed_maxs = results_3d['3d_maxs']
    threed_means = results_3d['3d_means']
    
    
    fig1, ax1 = plt.subplots(1, 3, figsize=(10, 3))
    fig2, ax2 = plt.subplots(1, 3, figsize=(10, 3))
    # fig1.suptitle("Summary Comparison of 0D and 3D Values at Relevant Points")

    s = 10

    # Reference 3D model Resultss
    ax1[1].scatter(range(len(threed_means)), threed_means, s = s, color = "black", marker='^', label = '3D')
    ax1[1].set_title("Mean",fontsize=fs)
    ax1[1].set_ylabel("Pressure [mmHg]",fontsize=fs)
    ax1[1].set_xlabel("Points",fontsize=fs)
    ax1[1].tick_params(axis='both', which='major', labelsize=fs)
    
    ax1[0].scatter(range(len(threed_maxs)), threed_maxs, s = s, color = "black", marker='^',label = '3D')
    ax1[0].set_title("Systolic",fontsize=fs)
    ax1[0].set_ylabel("Pressure [mmHg]",fontsize=fs)
    ax1[0].set_xlabel("Points",fontsize=fs)
    ax1[0].tick_params(axis='both', which='major', labelsize=fs)
    
    ax1[2].scatter(range(len(threed_mins)), threed_mins,s = s, color = "black", marker='^', label = '3D')
    ax1[2].set_title("Diastolic",fontsize=fs) 
    ax1[2].set_ylabel("Pressure [mmHg]",fontsize=fs)
    ax1[2].set_xlabel("Points",fontsize=fs)
    ax1[2].tick_params(axis='both', which='major', labelsize=fs)    

    
    # fig2.suptitle("3D vs 0D Pressures at Relevant Points.")

    ax2[0].plot([min(threed_maxs), max(threed_maxs)], [min(threed_maxs), max(threed_maxs)], 'k--', lw=0.8)
    ax2[0].set_title("Systolic",fontsize=fs)
    ax2[0].tick_params(axis='both', which='major', labelsize=fs)    

    ax2[1].plot([min(threed_means), max(threed_means)], [min(threed_means), max(threed_means)], 'k--', lw=0.8)
    ax2[1].set_title("Mean",fontsize=fs)
    ax2[1].tick_params(axis='both', which='major', labelsize=fs)    

    ax2[2].plot([min(threed_mins), max(threed_mins)], [min(threed_mins), max(threed_mins)], 'k--', lw=0.8)
    ax2[2].set_title("Diastolic",fontsize=fs)
    ax2[2].tick_params(axis='both', which='major', labelsize=fs)    

    for i in range(3):
        ax2[i].set_ylabel("0D Pressure [mmHg]",fontsize=fs)
        ax2[i].set_xlabel("3D Pressure [mmHg]",fontsize=fs)
    

    zerod = np.load(zerod_data, allow_pickle=True).item()

    markers = ['o','>','D']
    colors = ['b','r', 'm']
    
    for i, (name, results_0d) in enumerate(zerod.items()):
        zerod_mins = results_0d['0d_mins']
        zerod_maxs = results_0d['0d_maxs']
        zerod_means = results_0d['0d_means']
        ## Summary Statistics
        
        ## means
        ax1[1].scatter(range(len(zerod_means)), zerod_means,s = s, label = name, marker=markers[i], alpha = .7)        
        ## systolic
        ax1[0].scatter(range(len(zerod_maxs)), zerod_maxs, s = s,label = name, marker=markers[i], alpha = .7)        
        # diastolic
        ax1[2].scatter(range(len(zerod_mins)), zerod_mins,s = s, label = name, marker=markers[i], alpha = .7)
        
        ## plot 3D on x axis, and 0D on y axis
        s2=20
        ax2[0].scatter(threed_maxs, zerod_maxs,s = s2, label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)
        ax2[1].scatter(threed_means, zerod_means,s = s2, label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)
        ax2[2].scatter(threed_mins, zerod_mins, s = s2,label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)

    ax1[2].set_ylim([18.0,20.0])
    ax1[2].legend(fontsize=fs-2)
    ax2[2].legend(fontsize=fs-2)

    # handles, labels = ax1[0].get_legend_handles_labels()
    # fig1.legend(handles, labels, loc='upper right',fontsize=fs-2)
    # handles, labels = ax2[0].get_legend_handles_labels()
    # fig2.legend(handles, labels, loc='upper right',fontsize=fs-2)
    
    return fig1, fig2

if __name__ == '__main__':

    fig1, fig2 = plot_valid2('data/3d_data_f6.npy', 'data/0d_data_f6.npy')
    fig3, fig4 = plot_valid2('data/3d_data_f7.npy', 'data/0d_data_f7.npy')
    
    fig1.tight_layout()
    fig1.savefig("plots/f6.1.pdf")
    fig2.tight_layout()
    fig2.savefig("plots/f6.2.pdf")
    fig3.tight_layout()
    fig3.savefig("plots/f7.1.pdf")
    fig4.tight_layout()
    fig4.savefig("plots/f7.2.pdf")

