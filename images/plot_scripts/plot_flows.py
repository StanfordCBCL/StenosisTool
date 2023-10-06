from svinterface.core.bc import Inflow
from svinterface.plotting.params import set_params
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    set_params(use_latex=False, small_ticks=True)
    
    i1 = Inflow.from_file("data/diseased/AS1_SU0308_stent/flow_files/inflow_1D_orig.flow")
    i2 = Inflow.from_file("data/diseased/AS1_SU0308_stent/flow_files/inflow_1D.flow")
    
    ## save inflow graph
    fig,ax = plt.subplots(1,1, figsize = (8, 6))
    ax.plot(i1.t, i1.Q, 'r', label = 'Original Flow')
    ax.plot(i2.t, i2.Q, 'b', label = 'Corrected Flow')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Flow (ml/s)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend()
    
    plt.show()