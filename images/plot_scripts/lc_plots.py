import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.ticker import MaxNLocator


from pathlib import Path
from svinterface.core.polydata import Centerlines
from svinterface.plotting import params 
from svinterface.core.zerod.lpn import LPN


def load_data(c_3d: Centerlines, c_1d: Centerlines, lpn: LPN):
    """ load LC data"""
    
    # use valid array
    caps = c_3d.get_pointdata_array("Caps_0D")
    juncs = c_3d.get_pointdata_array("Junctions_0D") 
    vess = c_3d.get_pointdata_array("Vessels_0D") 
    #! pull out which outlet it actually is
    valid = np.array(sorted(list(set([0] + list(np.where(caps != -1)[0]) + list(np.where(juncs != -1)[0]) + list(np.where(vess != -1)[0])))))
    
    #! Remove later
    # gid_map = {}
    # tree = lpn.get_tree()
    # lpn.det_lpa_rpa(tree)
    # for node in lpn.tree_bfs_iterator(tree, allow='branch'):
    #     for v in node.vessel_info:
    #         for gid in v['gid']:
    #             gid_map[gid] = v['side']
    
    results_3d = {}
    # iterate through each valid point
    for oidx, point_id in enumerate(valid):
        time = []
        pressure = []
        flow = []
        #! use the fact that they should be in order already
        for arr_name in c_3d.get_pointdata_arraynames():
            if arr_name.startswith("pressure_"):
                time.append(float(arr_name.split('_')[1]))
                pressure.append(c_3d.polydata.GetPointData().GetArray(arr_name).GetValue(point_id))
        
        results_3d[oidx] = {'time': time,
                            'pressure': pressure,
                            'flow': flow,
                            'point_id': point_id}
    #     print(oidx, point_id, gid_map[point_id])
        
    # exit()
        
    results_1d = {}
    # iterate through each outlet
    for oidx, point_id in enumerate(valid):
        time = []
        pressure = []
        flow = []
        #! use the fact that they should be in order already
        for arr_name in c_1d.get_pointdata_arraynames():
            if arr_name.startswith("pressure_"):
                time.append(float(arr_name.split('_')[1]))
                pressure.append(c_1d.polydata.GetPointData().GetArray(arr_name).GetValue(point_id))
        
        results_1d[oidx] = {'time': time,
                            'pressure': pressure,
                            'flow': flow,
                            'point_id': point_id}

    zerod_means = []
    threed_means = []
    for i in range(len(valid)):
        zerod_means.append(np.trapz(results_1d[i]['pressure'], results_1d[i]['time']) / (results_1d[i]['time'][-1] - results_1d[i]['time'][0]))
        threed_means.append(np.trapz(results_3d[i]['pressure'], results_3d[i]['time']) / (results_3d[i]['time'][-1] - results_3d[i]['time'][0]))

    ## systolic
    zerod_maxs = []
    threed_maxs = []
    for i in range(len(valid)):
        zerod_maxs.append(np.array(results_1d[i]['pressure']).max())
        threed_maxs.append(np.array(results_3d[i]['pressure']).max())
   
    # diastolic
    zerod_mins = []
    threed_mins = []
    for i in range(len(valid)):
        zerod_mins.append(np.array(results_1d[i]['pressure']).min())
        threed_mins.append(np.array(results_3d[i]['pressure']).min())
        
    threed_results = {'3d_mins': threed_mins,
                      '3d_means': threed_means,
                      '3d_maxs': threed_maxs}
    
    zerod_results = {'0d_mins': zerod_mins,
                     '0d_means': zerod_means,
                     '0d_maxs': zerod_maxs}
    return threed_results, zerod_results
    

def rmse(x, xhat):
    
    return (((np.array(x) - np.array(xhat))**2).sum() / len(x))**(1/2)

def plot_valid(threed_data, zerod_data, names):
    
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

    ax2[0].plot([min(threed_maxs), max(threed_maxs)], [min(threed_maxs), max(threed_maxs)], 'k--', lw=0.8, label= 'Expected')
    ax2[0].set_title("Systolic",fontsize=fs)
    ax2[0].tick_params(axis='both', which='major', labelsize=fs)    

    ax2[1].plot([min(threed_means), max(threed_means)], [min(threed_means), max(threed_means)], 'k--', lw=0.8,label= 'Expected')
    ax2[1].set_title("Mean",fontsize=fs)
    ax2[1].tick_params(axis='both', which='major', labelsize=fs)    

    ax2[2].plot([min(threed_mins), max(threed_mins)], [min(threed_mins), max(threed_mins)], 'k--', lw=0.8,label= 'Expected')
    ax2[2].set_title("Diastolic",fontsize=fs)
    ax2[2].tick_params(axis='both', which='major', labelsize=fs)    

    for i in range(3):
        ax2[i].set_ylabel("0D Pressure [mmHg]",fontsize=fs)
        ax2[i].set_xlabel("3D Pressure [mmHg]",fontsize=fs)
    

    

    markers = ['o','>','D']
    colors = ['b','r', 'g']
    
    for i, (name, zerod_f) in enumerate(zip(names, zerod_data)):
        results_0d = np.load(zerod_f, allow_pickle=True).item()
        zerod_mins = results_0d['0d_mins']
        zerod_maxs = results_0d['0d_maxs']
        zerod_means = results_0d['0d_means']
        ## Summary Statistics
        
        ## means
        ax1[1].scatter(range(len(zerod_means)), zerod_means,s = s, color = colors[i], label = name, marker=markers[i], alpha = .7)        
        ## systolic
        ax1[0].scatter(range(len(zerod_maxs)), zerod_maxs, s = s,  color = colors[i], label = name, marker=markers[i], alpha = .7)        
        # diastolic
        ax1[2].scatter(range(len(zerod_mins)), zerod_mins,s = s, color = colors[i], label = name, marker=markers[i], alpha = .7)
        
        ## plot 3D on x axis, and 0D on y axis
        s2=20
        print(f"For {name}")
        ax2[0].scatter(threed_maxs, zerod_maxs,s = s2, label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)
        print(f"Corr Co (Systolic): {np.corrcoef(np.array(zerod_maxs), np.array(threed_maxs))[0,1]}")
        ax2[1].scatter(threed_means, zerod_means,s = s2, label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)
        print(f"Corr Co (Mean): {np.corrcoef(np.array(zerod_means), np.array(threed_means))[0,1]}")
        ax2[2].scatter(threed_mins, zerod_mins, s = s2,label = name, marker=markers[i], color=colors[i], alpha = .7, lw=0)
        print(f"Corr Co (Diastolic): {np.corrcoef(np.array(zerod_mins), np.array(threed_mins))[0,1]}")

    ax1[2].set_ylim([17.8,20.2])
    ax1[2].yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1[2].legend(fontsize=fs-2)
    ax2[2].legend(fontsize=fs-2)

    # handles, labels = ax1[0].get_legend_handles_labels()
    # fig1.legend(handles, labels, loc='upper right',fontsize=fs-2)
    # handles, labels = ax2[0].get_legend_handles_labels()
    # fig2.legend(handles, labels, loc='upper right',fontsize=fs-2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2

def plot_valid2(threed_data, zerod_data, names):
    markers = ['o','>','D']
    colors = ['r','b', 'g']
    fig1, ax1 = plt.subplots(3, 3, figsize=(10, 6), sharex=True)
    for idx, (data_3d, data_0d, name) in enumerate(list(zip(threed_data, zerod_data, names))):
        results_3d = np.load(data_3d, allow_pickle=True).item()
        threed_mins = results_3d['3d_mins']
        threed_maxs = results_3d['3d_maxs']
        threed_means = results_3d['3d_means']
        
        
        
        s = 10
        # Reference 3D model Resultss
        ax1[idx][1].scatter(range(len(threed_means)), threed_means, s = s, color = "black", marker='^', label = '3D ' + name )
        ax1[idx][1].tick_params(axis='both', which='major', labelsize=fs)
    
        ax1[idx][0].scatter(range(len(threed_maxs)), threed_maxs, s = s, color = "black", marker='^',label = '3D ' + name)
        ax1[idx][0].set_ylabel("Pressure [mmHg]",fontsize=fs)
        ax1[idx][0].tick_params(axis='both', which='major', labelsize=fs)
        
        ax1[idx][2].scatter(range(len(threed_mins)), threed_mins,s = s, color = "black", marker='^', label = '3D ' + name)
        ax1[idx][2].tick_params(axis='both', which='major', labelsize=fs)    

        if idx == 0:
            ax1[idx][1].set_title("Mean",fontsize=fs)
            ax1[idx][0].set_title("Systolic",fontsize=fs)
            ax1[idx][2].set_title("Diastolic",fontsize=fs) 
        
        if idx == 2:
            for i in range(3):
                ax1[idx][i].set_xlabel("Points",fontsize=fs)
        
        results_0d = np.load(data_0d, allow_pickle=True).item()
        zerod_mins = results_0d['0d_mins']
        zerod_maxs = results_0d['0d_maxs']
        zerod_means = results_0d['0d_means']
        ## Summary Statistics
        
        ## means
        ax1[idx][1].scatter(range(len(zerod_means)), zerod_means,s = s, label ='0D ' + name, marker=markers[0],color=colors[idx], alpha = .7)        
        ## systolic
        ax1[idx][0].scatter(range(len(zerod_maxs)), zerod_maxs, s = s,label ='0D ' + name, marker=markers[0],color=colors[idx], alpha = .7)        
        # diastolic
        ax1[idx][2].scatter(range(len(zerod_mins)), zerod_mins,s = s, label = '0D ' + name, marker=markers[0], color=colors[idx],alpha = .7)
            

        ax1[idx][2].set_ylim([17.8,20.2])
        ax1[idx][2].yaxis.set_major_locator(MaxNLocator(integer=True))
        ax1[idx][2].legend(fontsize=fs-2)

 
    fig1.tight_layout()
    return fig1

if __name__ == '__main__':
    
    # Parameter for loading from npy if already saved
    # quick_load = False
    quick_load = False
    
    fs=12
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.rc('text', usetex=True)
    params.set_params(small_ticks=True, use_latex=True, plw = 0.01)

    dis_filenames = ['no_lc', 'original_lc', 'split_lc', None, 'all_lc', 'iter_lc']
    
    if not quick_load:
        
        # load the data from source & save it
        c_3d_dis = Centerlines.load_centerlines("data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/3D_DIR/prestent/AS1_SU0308_3D_centerlines.formatted.vtp")
        lpn_dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/LPN_DIR')
        for sim, name in enumerate(dis_filenames):
            if name is None:
                continue
            print(f"Loading {name}...", end = '', flush = True)
            new_dir = lpn_dir / f'AS1_SU0308.sim.{sim}'
            
            c1d = Centerlines.load_centerlines(str(new_dir / 'centerline_projection.vtp'))
            #lpn = LPN.from_file(str(new_dir / "AS1_SU0308.in"))
            threed_results, zerod_results = load_data(c_3d=c_3d_dis, c_1d=c1d)#, lpn=lpn)
            if sim == 0:
                np.save('images/plot_data/threed.npy', threed_results, allow_pickle=True)
            np.save(f'images/plot_data/{name}.npy', zerod_results, allow_pickle=True)
            print("Done")
        
        # load the repaired data
        param_dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/parameterization')
        rep_3d_dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/3D_DIR')
        for sim, (name, threed_name) in enumerate(zip(['LPA_all', 'LPA_limited', 'RPA_all', 'RPA_limited', 'RPA_2_all', 'RPA_2_limited'],
                                       ['LPA_stent/AS1_SU0308_3D_LPA_stented_centerlines.mapped.formatted.vtp','LPA_stent/AS1_SU0308_3D_LPA_stented_centerlines.mapped.formatted.vtp', 'RPA_stent/AS1_SU0308_3D_RPA_stented_centerlines.mapped.formatted.vtp','RPA_stent/AS1_SU0308_3D_RPA_stented_centerlines.mapped.formatted.vtp', 'RPA_2_stent/AS1_SU0308_3D_RPA_stented_2_centerlines.mapped.formatted.vtp', 'RPA_2_stent/AS1_SU0308_3D_RPA_stented_2_centerlines.mapped.formatted.vtp'])):
            print(f"Loading {name}...", end = '', flush = True)
            rep_3d_file = rep_3d_dir / threed_name.split('/')[0] / threed_name.split('/')[1]
            rep_c3d = Centerlines.load_centerlines(str(rep_3d_file))
            rep_c1d = Centerlines.load_centerlines(str(param_dir / name / 'centerline_projection.vtp'))
            
            threed_results, zerod_results = load_data(c_3d=rep_c3d, c_1d=rep_c1d)
            np.save(f'images/plot_data/{name}_3d.npy', threed_results, allow_pickle=True)
            np.save(f'images/plot_data/{name}.npy', zerod_results, allow_pickle=True)
            print("Done")
        
    
    # Fig 7, 8
    dis_3d_file = 'images/plot_data/threed.npy'
    
    dis_0d_files = ['images/plot_data/' + file + '.npy' for file in  ['no_lc', 'original_lc']]
    fig1, fig2 = plot_valid(dis_3d_file, dis_0d_files, names = ['0D No Correction', '0D One-Step Correction'])
    
    dis_0d_files = ['images/plot_data/' + file + '.npy' for file in  ['split_lc', 'all_lc', 'iter_lc']]
    fig3, fig4 = plot_valid(dis_3d_file, dis_0d_files,  names = ['0D Subdomains', '0D S + Vessels', '0D S + V + Iterative'])
    exit()
    fig1.savefig("images/paper/04_zerod/a_lc_nlc.pdf")
    fig2.savefig("images/paper/04_zerod/b_lc_nlc.pdf")
    fig3.savefig("images/paper/04_zerod/a_lc_ext.pdf")
    fig4.savefig("images/paper/04_zerod/b_lc_ext.pdf")
    
    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)
    plt.close(fig4)
    
    
    dis_3d_files = ['images/plot_data/' + file + '.npy' for file in ['LPA_limited_3d', 'RPA_limited_3d', 'RPA_2_limited_3d']]
    dis_0d_files = ['images/plot_data/' + file + '.npy' for file in  ['LPA_limited', 'RPA_limited', 'RPA_2_limited']]
    fig5 = plot_valid2(dis_3d_files, dis_0d_files, names = ['LPA Proximal', 'RPA Proximal', 'RPA Distal'])
    fig5.savefig("images/paper/05_param/repaired_lc.pdf")