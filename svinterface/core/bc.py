import os
from collections import OrderedDict
import re
import numpy as np
from scipy.interpolate import interp1d

class RCR(object):
    '''The RCR class is used to manipulate 0D RCR simulation boundary conditions with limited in-build 3D-0D RCR convertability.
    
    
    Attributes:
        bc_list (list[dict]): The list of boundary conditions.
    '''

    # BC types.
    BC_TYPE_RCR = "RCR"
    BC_TYPE_RESISTANCE = "Resistance"
    BC_TYPE_PRESCRIBED_VELOCITIES = "Prescribed Velocities"

    # File names storing BC values for each BC type.
    RCR_FILE_NAME = "rcrt.dat"
    RESISTANCE_FILE_NAME = "resistance.dat"

    def __init__(self):
        self.bc_list = OrderedDict()
        self.valid_names = True

    def add_rcr(self, face_name, Rp, C, Rd, Pd):
        self.bc_list[face_name] = { 'type': self.BC_TYPE_RCR, 'faceID': face_name, 'Rp':Rp, 'C':C, 'Rd':Rd, 'Pd': Pd}

    def write_rcrt_file(self, dirpath=None, as_3d = False):
        '''Write RCR boundary conditions to a file.
        sort_for_3d can be used prior to writing to ensure ordering.
        '''
        # check there are BC
        num_bcs = sum([bc['type'] == self.BC_TYPE_RCR for bc in self.bc_list])
        if num_bcs == 0:
            return
        
        # check that filenames are valid for 0D
        if not as_3d and not self.valid_names:
            raise KeyError("Cannot write as 0D file without valid mapping of faceID's to rcr's")

        bc_path = dirpath

        newline = os.linesep 
        with open(bc_path + os.sep + self.RCR_FILE_NAME, "w") as rcr_file:
            rcr_file.write('2' + newline)
            for bc in self.bc_list:
                if bc['type'] != self.BC_TYPE_RCR:
                    continue
                rcr_file.write('2' + newline)
                if not as_3d:
                    rcr_file.write(str(bc['faceID']) + newline) 
                for pname in ['Rp', 'C', 'Rd']:
                    rcr_file.write(str(bc[pname]) + newline) 
                pressure = str(bc['Pd'])
                rcr_file.write('0.0 ' + pressure + newline) 
                rcr_file.write('1.0 ' + pressure + newline) 

    def read_rcrt_file(self, rcrt_file, as_3d = False, solver_inp: str = None, svpre: str = None):
        ''' Read RCR BC from file. Assigned arbitrary counter for 3D files
        '''
        # check if mapped
        mapped = not as_3d
        if as_3d and solver_inp and svpre:
            _, rcrt_list = self._parse_svsolver(solver_inp)
            id_to_name = self._parse_svpre(svpre)
            mapped = True
        
        # read file
        with open(rcrt_file, 'r') as rfile:
            keyword = rfile.readline()
            counter = 0
            while True:
                tmp = rfile.readline()
                
                if tmp == keyword:
                    if as_3d:
                        face_name = id_to_name[rcrt_list[counter]] if mapped else counter
                    else:
                        face_name = rfile.readline().rstrip()
                    Rp = float(rfile.readline())
                    C = float(rfile.readline())
                    Rd = float(rfile.readline())
                    p0 = float(rfile.readline().strip().split()[1])
                    p1 = float(rfile.readline().strip().split()[1])
                    assert p0 == p1, 'Cannot handle time-dependent reference pressure'
                    Pd = (float(p1))
                    
                    self.add_rcr(face_name = face_name, Rp = Rp, C = C, Rd = Rd, Pd = Pd)
                if len(tmp) == 0:
                    break
                counter += 1
        
        # set validity
        self.valid_names = mapped
        return

    def sort_for_3d(self, solver_inp: str, svpre: str):
        """sorts boundary conditions for 3D use from 0D

        Args:
            solver_inp (str): path to a 3D simulation's solver.inp file
            svpre (str): path to a 3D simulation's svpre file
            
        """
        
        # parse files
        rcrt_num, rcrt_list = self._parse_svsolver(solver_inp)
        id_to_name = self._parse_svpre(svpre)
        
        assert rcrt_num == len(self.bc_list), "Number of RCR Surfaces doesn't match number of BC's listed. "
        
        # reorder bc_list
        bc_list = OrderedDict()
        
        if self.valid_names:
            # if valid_names, then sort according to name (sorting 0D to 3D)
            for idx in rcrt_list:
                bc_list[id_to_name[idx]] = self.bc_list[id_to_name[idx]]
            self.bc_list = bc_list
        else:
            raise Exception("Sort is not guarenteed since no faces were mapped.")
        
    def _parse_svpre(self, svpre: str):
        """Parses an svPre file to get id of 3D simulation"""
        id_to_name = {}
        with open(svpre, 'r') as sfile:
            for line in sfile:
                match= re.search("set_surface_id_vtp mesh-complete/mesh-surfaces/(.*).vtp ([0-9]*)", line)
                if match:
                        id_to_name[int(match.group(2))] = match.group(1)
                    
        return id_to_name
    
    def _parse_svsolver(self, inp: str):
        ''' parses an solver.inp file to extract the order of RCRTs in the rcrt file
        
        Returns:
            (int, list): (number of RCRs)
        '''
        with open(inp, 'r') as ifile:
            f = ''.join(ifile.readlines())
            rcr_num = int(re.search("^Number of RCR Surfaces: ([0-9]*)$", f, re.MULTILINE).group(1))
            rcr_list = (re.search("^List of RCR Surfaces: (.*)$", f, re.MULTILINE).group(1).split(' '))
            rcr_list = [int(x) for x in rcr_list]
        return rcr_num, rcr_list
    


class Inflow():
    ''' Handles inflow and inflow files
    '''
    def __init__(self, inflow_arr, inverse = False, smooth = True, n_points = 1000):
        ''' inflow arr must be inform np.array[(t, Q), (t, Q), ...] 
        '''
        self.inflow = inflow_arr
        
        self.correct_flow()
        
        if inverse:
            self.inverse_flow()
        if smooth:
            self.smooth_flow(n_points)
        
        self.update()
       
    
    def update(self):
        """Updates internal attributes when operations change values
        """
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
        self.update()

    def correct_flow(self):
        ''' Checks that the first and last flow values are equivalent. If not, append a last term that is equivalent to the first term
        '''
        if self.inflow[0, 1] - self.inflow[-1, 1] != 0:
            time_diff = self.inflow[1, 0] - self.inflow[0,0]
            self.inflow = np.append(self.inflow, np.array([[time_diff + self.inflow[-1, 0], self.inflow[0, 1]]]), axis = 0)
        self.update()
            
    def smooth_flow(self, n_points):
        ''' smooth flow using a cubic spline 
        '''
        f = interp1d(self.inflow[:, 0], self.inflow[:, 1], kind = 'cubic')
        x = np.linspace(self.inflow[0, 0], self.inflow[-1, 0], n_points)
        y = f(x)
        self.inflow = np.array(list(zip(x, y)))
        self.update()
          