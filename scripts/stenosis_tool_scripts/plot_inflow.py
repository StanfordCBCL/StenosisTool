import matplotlib.pyplot as plt
from functions.flow import Inflow
import sys
import os

if __name__ == '__main__':
    
    iflow = Inflow(sys.argv[1], smooth = False)
    
    ## save inflow graph
    fig,ax = plt.subplots(1,1 )
    ax.plot(iflow.t, iflow.Q)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('flow (ml/s)')
    ax.set_title('Inflow')
    fig.savefig(os.path.join(os.path.dirname(sys.argv[1]), os.path.splitext(os.path.basename(sys.argv[1]))[0] + '.png'))