import os

class BoundaryConditions(object):
        '''The BoundaryConditions class is used to set 0D simulation boundary conditions.
        Attributes:
            bc_list (list[dict]): The list of boundary conditions.
            bc_path (str): The path to the boundary conditions files.
        '''

        # BC types.
        BC_TYPE_RCR = "RCR"
        BC_TYPE_RESISTANCE = "Resistance"
        BC_TYPE_PRESCRIBED_VELOCITIES = "Prescribed Velocities"

        # File names storing BC values for each BC type.
        RCR_FILE_NAME = "rcrt.dat"
        RESISTANCE_FILE_NAME = "resistance.dat"

        def __init__(self):
            self.bc_list = []

        def add_resistance(self, face_name, resistance):
            self.bc_list.append( { 'type': self.BC_TYPE_RESISTANCE, 'faceID': face_name, 'resistance':resistance})

        def add_rcr(self, face_name, Rp, C, Rd, Pd):
            self.bc_list.append( { 'type': self.BC_TYPE_RCR, 'faceID': face_name, 'Rp':Rp, 'C':C, 'Rd':Rd, 'Pd': Pd})

        def add_velocities(self, face_name, file_name):
            self.bc_list.append( { 'type': self.BC_TYPE_PRESCRIBED_VELOCITIES, 'faceID': face_name, 'file_name': file_name} )


        def write_files(self, path=None):
            '''Write boundary conditions to files for each specific type.
            '''
            self.write_rcrt_file(path)
            self.write_resistance_file(path)

        def write_resistance_file(self, path=None):
            '''Write RESISTANCE boundary conditions to a file.
            '''
            num_bcs = sum([bc['type'] == self.BC_TYPE_RESISTANCE for bc in self.bc_list])
            if num_bcs == 0:
                return

            if path == None:
                bc_path = self.bc_path
            else:
                bc_path = path

            newline = os.linesep 
            with open(bc_path + os.sep + self.RESISTANCE_FILE_NAME, "w") as res_file:
                for bc in self.bc_list:
                    if bc['type'] != self.BC_TYPE_RESISTANCE:
                        continue
                    res_file.write(bc['faceID'] + ' ' + str(bc['resistance']) + newline) 


        def write_rcrt_file(self, path=None, three_d = False):
            '''Write RCR boundary conditions to a file.
            '''
            num_bcs = sum([bc['type'] == self.BC_TYPE_RCR for bc in self.bc_list])
            if num_bcs == 0:
                return

            if path == None:
                bc_path = self.bc_path
            else:
                bc_path = path

            newline = os.linesep 
            with open(bc_path + os.sep + self.RCR_FILE_NAME, "w") as rcr_file:
                rcr_file.write('2' + newline)
                for bc in self.bc_list:
                    if bc['type'] != self.BC_TYPE_RCR:
                        continue
                    rcr_file.write('2' + newline)
                    if not three_d:
                        rcr_file.write(str(bc['faceID']) + newline) 
                    for pname in ['Rp', 'C', 'Rd']:
                        rcr_file.write(str(bc[pname]) + newline) 
                    pressure = str(bc['Pd'])
                    rcr_file.write('0.0 ' + pressure + newline) 
                    rcr_file.write('1.0 ' + pressure + newline) 

        def read_rcrt_file(self, rcrt_file, three_d = False):
            ''' Read rcr BC from file
            '''
            with open(rcrt_file, 'r') as rfile:
                keyword = rfile.readline()
                while True:
                    tmp = rfile.readline()
                    
                    if tmp == keyword:
                        if three_d:
                            face_name = ''
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
            return
                    
        def get_bc_map(self):
            bc_map = {}
            for bc in self.bc_list:
                bc_map[bc['faceID']] = bc
            return bc_map
            
