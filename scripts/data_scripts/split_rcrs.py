from src.misc import create_parser
from src.data_org import DataPath

import numpy as np
import os

###########
# Classes #
###########

class BoundaryConditions(object):
        '''The BoundaryConditions class is used to set 1D simulation boundary conditions.
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

        def write_rcrt_file(self, path=None):
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
                    for pname in ['faceID', 'Rp', 'C', 'Rd']:
                        rcr_file.write(str(bc[pname]) + newline) 
                    pressure = str(bc['Pd'])
                    rcr_file.write('0.0 ' + pressure + newline) 
                    rcr_file.write('1.0 ' + pressure + newline) 


#############
# Split RCR #
#############

def load_area_file(area_filepath):
    with open(area_filepath, 'r') as afile:
        areas = {}
        afile.readline() # ignore first comment line
        for line in afile:
            line = line.rstrip().split()
            areas[line[0]] = float(line[1])
    return areas

def total_area(areas):
    return sum(list(areas.values()))

def validate_caps(areas):
    ''' confirm that all caps have either lpa or rpa in them'''
    names = list(areas.keys())
    for name in names:
        if 'lpa' not in name.lower() and 'rpa' not in name.lower():
            raise ValueError('Unable to identify RPA vs. LPA caps: please rename caps to contain lpa or rpa')
    return 

def split_rpa_lpa(areas):
    ''' splits areas between lpa and rpa areas '''
    lpa_areas = {}
    rpa_areas = {}
    
    validate_caps(areas)
    
    for name, area in areas.items():
        if 'lpa' in name.lower():
            lpa_areas[name] = area
        elif 'rpa' in name.lower():
            rpa_areas[name] = area
    return lpa_areas, rpa_areas


def split_bc(areas,  x):
    
    def Rpi(Ai, A, Rp):
        return (A / Ai) * Rp
    
    def Rdi(Ai, A, Rd):
        return (A / Ai) * Rd
    
    def Ci(Ai, A, C):
        return (Ai/A) * C
    
    lpa_areas, rpa_areas = split_rpa_lpa(areas)
    
    lpa_A = total_area(lpa_areas)
    rpa_A = total_area(rpa_areas)
    
    rcrs = {}
    
    
    for name, area in lpa_areas.items():
        rcrs[name] = {'Rp': Rpi(area, lpa_A, x[0]),
                     'C': Ci(area, lpa_A, x[1]),
                     'Rd': Rdi(area, lpa_A, x[2]) }
    
    for name, area in rpa_areas.items():
        rcrs[name] = {'Rp': Rpi(area, rpa_A, x[3]),
                     'C': Ci(area, rpa_A, x[4]),
                     'Rd': Rdi(area, rpa_A, x[5]) }
    
    return rcrs

def calc_rcrs(results, inlet, area_file, out_dir):
    ''' Splits RCRS '''
    
    areas = load_area_file(area_file)
    
    del areas[inlet]
    
    rcrs = split_bc(areas, results['x'])
    
    bcs = BoundaryConditions()
    
    for name, vals in rcrs.items():
        bcs.add_rcr(face_name=name, Rp = vals['Rp'], C = vals['C'], Rd = vals['Rd'], Pd = results['Pd'])
    
    bcs.write_rcrt_file(out_dir)
    
    return rcrs


'''
def convert_old_rcrt(inlet, mdl_cvpre_file, old_rcrt_file, solver3d, out_dir, ingrid = False):
    
    # map id to name
    if not ingrid:
        face_mappings = parse_mdl(mdl_cvpre_file, reverse = True)
    else:
        face_mappings = ingrid_rcrt_map(mdl_cvpre_file)
        
    # get used rcrt values
    with open(solver3d, 'r') as sfile:
        for line in sfile:
            if re.search('List of RCR Surfaces:', line):
                used_values = line.split(':')[1].split()
                break
    
    print(used_values, face_mappings)
        
    with open(old_rcrt_file, 'r') as rfile:
        keyword = rfile.readline()
        cur_cap_idx = 0
        bcs = BoundaryConditions()
        while True:
            #print(cur_cap_idx, face_mappings[ids[cur_cap_idx]])
            
            
            tmp = rfile.readline()
            if tmp == keyword:
                face_name = face_mappings[int(used_values[cur_cap_idx])]
                Rp = float(rfile.readline())
                C = float(rfile.readline())
                Rd = float(rfile.readline())
                p0 = float(rfile.readline().strip().split()[1])
                p1 = float(rfile.readline().strip().split()[1])
                assert p0 == p1, 'Cannot handle time-dependent reference pressure'
                Pd = (float(p1))
                
                # add rcr bc
                bcs.add_rcr(face_name = face_name, Rp = Rp, C = C, Rd = Rd, Pd = Pd)
            
            if cur_cap_idx > len(used_values):
                raise RuntimeError('More rcrt BCs than possible caps')
            
            cur_cap_idx += 1
            if len(tmp) == 0:
                break
            
        if cur_cap_idx <= len(used_values):
            raise RuntimeError('Insuffienct rcrt BCs to match caps')
                
    bcs.write_rcrt_file(out_dir)


def ingrid_rcrt_map(cvpre_file):
    with open(cvpre_file, 'r') as cvprefile:
        while True:
            tmp = cvprefile.readline()
            if tmp == '# Assign IDs to the surfaces\n':
                break
        mapping = {}
        while True:
            tmp = cvprefile.readline().rstrip().split()
            if tmp == []:
                break
            cap_name = 'cap_' + os.path.splitext(os.path.basename(tmp[1]))[0]
            cap_id = int(tmp[2])
            mapping[cap_id]= cap_name
    print(mapping)
    return mapping
'''

#########
# Mains #
#########

def tool_main(args):
    # TODO: CREATE MODEL AND SPLIT
    raise NotImplementedError

def dev_main(args):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        if 'cap_info' in model.info['files']:
            area_file = model.info['files']['cap_info']
        else:
            print('Missing CapInfo for ' + model_name)
            continue
        results_file = os.path.join(os.path.dirname(model.tune_solver), model.info['metadata']['name'] + '_optimize_results.npy')
        results = np.load(results_file, allow_pickle=True ).item()
        
        calc_rcrs(results=results, 
                  inlet=model.info['model']['inlet'],
                  area_file=area_file,
                  out_dir=model.solver_dir)
            


if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'Splits optimized results in rcrt.dat')
    
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)