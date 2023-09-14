# File: evaluate_nn.py
# File Created: Tuesday, 15th August 2023 12:44:58 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:39:02 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Evaluates the Neural Network


import matplotlib.pyplot as plt
import torch
import tqdm
from pathlib import Path
import numpy as np

from svinterface.plotting.params import set_params
from svinterface.utils.misc import d2m

def plot_avp(y, savedir: Path):
    """ Plot Actual vs Predicted """
    # plot for all values at a point
    for i in tqdm.tqdm(range(0,len(y[0][0]),6), desc='Generating A vs. P plots'):
        fig, ax = plt.subplots(2, 3, figsize = (24, 16), )
        ax = ax.flatten()
        names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)', 'Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
        for j in range(6):
            if j < 3:
                ax[j].scatter(y[:,0,i+j],y[:,1,i+j])
            else:
                ax[j].scatter(d2m(y[:,0,i+j]),d2m(y[:,1,i+j]))
            ax[j].set_xlabel("Actual")
            ax[j].set_ylabel("Predicted")
            ax[j].set_title(names[j])
        fig.suptitle(f"Truth vs. Predicted at point {i // 6 +1}")
        fig.savefig(str(savedir / f'avp_point_{i//6+1}.png'))
        plt.close(fig)

def plot_overlap(x, y, savedir: Path):
    """ Plot Scatter plots"""

    indices = np.random.choice(list(range(len(x))), size=len(x)//4)
    
    # plot for all values at a point
    for i in tqdm.tqdm(range(0,len(y[0][0]),6), desc='Generating Scatter Plots'):
        fig, ax = plt.subplots(3, 2, figsize = (16, 24), )
        for row in range(3):
            for col, name  in enumerate(['Flow', 'Pressure']):
                if col == 0:
                    name2=name + ' (mL/s)'
                    units = lambda x:x
                else:
                    name2=name + ' (mmHg)'
                    units = d2m
                
                names = [f'Diastolic {name}', f'Mean {name}', f'Systolic {name}']
                colors = ['red', 'blue', 'green']
                for j in range(3):
                    ax[row][col].scatter(x[indices, row], units(y[indices, 0, i + j + 3*col]),c = colors[j], marker = '^', alpha=.5, label='Actual ' + names[j])
                    ax[row][col].scatter(x[indices, row], units(y[indices, 1, i + j + 3*col]),c = colors[j], alpha=.5, label='Predicted ' + names[j])
                
                ax[row][col].set_xlabel(f'Repair Param {row + 1}')
                ax[row][col].set_ylabel(name2)
                ax[row][col].set_title(f'Repair Param {row + 1} {name}')
        
    
        fig.suptitle(f"Scatter at point {i // 6 +1}")
        fig.savefig(str(savedir / f'scatter_point_{i//6+1}.png'))
        plt.close(fig)
    
def print_stats(ys):
    """ Print Model Statistics """
    residuals = ys[:,0] - ys[:,1]
    r = np.abs(residuals / ys[:,0])
    ABSRELE = r.mean().item()
    MAXRELE = r.max().item()
    sd = r.std().item()
    print(f"Mean Relative Error: {round(ABSRELE * 100, 4)}%" )
    print(f"Maximum Relative Error: {round(MAXRELE * 100, 4)}%")
    print(f"2 std Relative Error: (0%, {round(ABSRELE + 2 * sd * 100, 4)}%)")
    
    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Evaluates the neural network")
    
    dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR')
    
    # load 
    version = 'version_0'
    x = torch.load(dir / "training_results" /  "run1" / "lightning_logs" / version / "predict_input.pt")
    ys = torch.load(dir / "training_results" /  "run1" / "lightning_logs" / version / "predict_output.pt")
    
    plots_dir = dir / "training_results" /  "run1" / "lightning_logs" / version / 'plots'
    
    # set matplotlib params
    set_params()
    
    print_stats(ys)
    
    avp_dir = plots_dir / 'actual_vs_predicted'
    avp_dir.mkdir(parents=True, exist_ok=True)
    
    plot_avp(ys, savedir=avp_dir)
    
    overlap_dir = plots_dir / 'overlap'
    overlap_dir.mkdir(parents=True, exist_ok=True)
    plot_overlap(x, ys, savedir=overlap_dir)