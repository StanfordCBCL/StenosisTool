from core.bc import Inflow
import matplotlib.pyplot as plt

def plot_flow(self, inflow: Inflow,  save = False, output_file = None):
    ''' plot the flow
    '''
    
    ## save inflow graph
    fig,ax = plt.subplots(1,1, figsize = (8, 6))
    ax.plot(inflow.t, inflow.Q)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Flow (ml/s)')
    fig.suptitle("Inflow Waveform")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    if save:
        fig.savefig(output_file)
    else:
        plt.show()
        plt.close(fig)