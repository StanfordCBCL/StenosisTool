import numpy as np
from scipy.interpolate import interp1d

class Inflow():
    ''' Handles inflow and inflow files'''
    def __init__(self, inflow_file, inverse = False, smooth = True, n_points = 1000):
        self.inflow = np.loadtxt(inflow_file)
        
        self.check_valid_flow()
        
        if inverse:
            self.inverse_flow()
        if smooth:
            self.smooth_flow(n_points)
        
        self.t = self.inflow[:, 0]
        self.Q = self.inflow[:, 1]
        # cycle length
        self.tc = (self.t[-1] - self.t[0])
        
        self.mean_inflow = np.trapz(self.Q, self.t) / self.tc
        self.max_inflow = self.Q.max()
        self.min_inflow = self.Q.min()
    
    def write_flow(self, flowpath):
        np.savetxt(flowpath, self.inflow)
    
    def inverse_flow(self):
        '''inverse the pos and negative flow values'''
        self.inflow[:, 1] *= -1

    def check_valid_flow(self):
        if self.inflow[0, 1] - self.inflow[-1, 1] != 0:
            time_diff = self.inflow[1, 0] - self.inflow[0,0]
            self.inflow = np.append(self.inflow, np.array([[time_diff + self.inflow[-1, 0], self.inflow[0, 1]]]), axis = 0)
            
            
    def smooth_flow(self, n_points):

        f = interp1d(self.inflow[:, 0], self.inflow[:, 1], kind = 'cubic')
        x = np.linspace(self.inflow[0, 0], self.inflow[-1, 0], n_points)
        y = f(x)
        
        self.inflow = np.array(list(zip(x, y)))
          