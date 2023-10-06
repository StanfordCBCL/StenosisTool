from svinterface.core.bc import Inflow
import matplotlib.pyplot as plt

def plot_flow(inflow: Inflow,  save = False, output_file = None):
    ''' plot the flow
    '''
    
    ## save inflow graph
    fig,ax = plt.subplots(1,1, figsize = (4, 3))
    ax.plot(inflow.t, inflow.Q)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Flow (ml/s)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    if save:
        fig.savefig(output_file)
    return fig