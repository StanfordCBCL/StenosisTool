import numpy as np
from scipy.interpolate import interp1d
class Inflow0D():
    ''' Handles inflow and inflow files'''
    
    
    def __init__(self, inflow_arr, inverse = False, smooth = True, n_points = 1000):
        ''' inflow arr must be inform np.array[(t, Q), (t, Q), ...] 
        '''
        self.inflow = inflow_arr
        
        self.correct_flow()
        
        if inverse:
            self.inverse_flow()
        if smooth:
            self.smooth_flow(n_points)
            
        # compute values
        self.t = self.inflow[:, 0]
        self.Q = self.inflow[:, 1]
        self.tc = (self.t[-1] - self.t[0])
        
        self.mean_inflow = np.trapz(self.Q, self.t) / self.tc
        self.max_inflow = self.Q.max()
        self.max_inflow_t = self.t[np.where(self.Q == self.max_inflow)][0]
        self.min_inflow = self.Q.min()
        self.min_inflow_t = self.t[np.where(self.Q == self.min_inflow)][0]
    
    
    @classmethod
    def from_file(cls, inflow_file, inverse = False, smooth = True, n_points = 1000):
        ''' Constructs inflow from a file
        '''
        inflow_arr = np.loadtxt(inflow_file,)
        return cls(inflow_arr, inverse, smooth, n_points)

    def write_flow(self, flow_path):
        ''' Writes inflow to a file
        '''
        np.savetxt(flow_path, self.inflow)
    
    def inverse_flow(self):
        '''inverse the pos and negative flow values
        '''
        self.inflow[:, 1] *= -1

    def correct_flow(self):
        ''' Checks that the first and last flow values are equivalent. If not, append a last term that is equivalent to the first term
        '''
        if self.inflow[0, 1] - self.inflow[-1, 1] != 0:
            time_diff = self.inflow[1, 0] - self.inflow[0,0]
            self.inflow = np.append(self.inflow, np.array([[time_diff + self.inflow[-1, 0], self.inflow[0, 1]]]), axis = 0)
            
            
    def smooth_flow(self, n_points):
        ''' smooth flow using a cubic spline 
        '''

        f = interp1d(self.inflow[:, 0], self.inflow[:, 1], kind = 'cubic')
        x = np.linspace(self.inflow[0, 0], self.inflow[-1, 0], n_points)
        y = f(x)
        self.inflow = np.array(list(zip(x, y)))
          
    def plot_flow(self, output_file):
        ''' plot the flow
        '''
        
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            print(e, ': aborting plot')
            return
        
        ## save inflow graph
        fig,ax = plt.subplots(1,1 )
        ax.plot(self.t, self.Q)
        ax.set_xlabel('time (s)')
        ax.set_ylabel('flow (ml/s)')
        fig.savefig(output_file)