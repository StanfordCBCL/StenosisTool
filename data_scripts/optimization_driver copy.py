from src.data_org import DataPath
from src.flow import Inflow
from src.centerlines import Centerlines
from src.file_io import Solver0D, 
from src.misc import m2d, d2m, create_parser

import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
import svzerodsolver as zerod
import os
import matplotlib.pyplot as plt
import json
from tqdm import tqdm
from functools import partialmethod
import re
from argparse import Namespace


# disable tqdm
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

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


class TuneParams():
    
    def __init__(self):
        ## optimizer params defaults
        self.cardiac_cycles = 6
        self.num_timesteps_per_cycle = 1000
        self.rpa_flow_split = .55
        self.mPAP_meas = m2d(15) # mmHg -> cgs
        self.cap_wedge_pressure = m2d(7) # mmHg -> cgs
        self.viscosity = .04 # cgs
        self.density = 1.06 # cgs
        self.linear_ehr = 1.2e6 # dyne/cm^2
        
        
        ## for termination
        self.pat = 5
        self.pat_tol  = 1e-6
        self.tolerance = 0.01


##################
# Functions
###################

def convert_to_dict(opt_results: optimize.OptimizeResult):
    rez = {}
    for key, val in opt_results.items():
        rez[key] = val
    return rez


def write_0d_dict(params: TuneParams, inflow: Inflow, centerlines: Centerlines, units = 'mm') -> Solver0D:
    ''' sets up solver file skeleton'''
    solver_data = Solver0D()
    solver_data.set_up_new_solver()
    
    ## bc
    # inlet
    inflow_bc= {'bc_name': 'INFLOW',
              'bc_type': 'FLOW',
              'bc_values': {'Q': list(inflow.Q),
                            't': list(inflow.t)}}
    # outlets
    lpa_bc = {'bc_name': 'LPA_BC',
              'bc_type': 'RESISTANCE',
              "bc_values": {"Pd": params.cap_wedge_pressure,
                            "R": 0}}
    
    rpa_bc = {'bc_name': 'RPA_BC',
              'bc_type': 'RESISTANCE',
              "bc_values": {"Pd": params.cap_wedge_pressure,
                            "R": 0}}
    
    solver_data.bc = [inflow_bc, lpa_bc, rpa_bc]
    
    
    ## junction basic
    solver_data.junctions=[
        {
            "inlet_vessels": [
                0
            ],
            "junction_name": "J0",
            "junction_type": "BloodVesselJunction",
            "outlet_vessels": [
                1,
                3
            ]
        }, 
        {
            "inlet_vessels": [
                1
            ],
            "junction_name": "J1",
            "junction_type": "internal_junction",
            "outlet_vessels": [
                2
            ]
        },
        {
            "inlet_vessels": [
                3
            ],
            "junction_name": "J2",
            "junction_type": "internal_junction",
            "outlet_vessels": [
                4
            ]
        },
    ]

    ## simulation_parameters
    solver_data.simulation_params = {
        "number_of_cardiac_cycles": params.cardiac_cycles,
        "number_of_time_pts_per_cardiac_cycle": params.num_timesteps_per_cycle,
        "density": params.density,
        "viscosity": params.viscosity,
    }
    
    ## vessels
    
    # Compute MPA
    if units == 'cm':
        coeff = 1
    elif units == 'mm': # convert to cgs
        coeff = .1
    else:
        raise ValueError('units much be cm or mm')
    paths = centerlines.get_pointdata(centerlines.PointDataFields.PATH) * coeff
    areas = centerlines.get_pointdata(centerlines.PointDataFields.AREA) * coeff**2
    branches = centerlines.get_pointdata(centerlines.PointDataFields.BRANCHID)
    
    points_id = np.where(branches == 0)
    
    length = paths[points_id][-1]
    # as opposed to only inlet and outlet, we use average of the entire section
    avg_radius = np.sqrt(areas[points_id]/np.pi).mean()
    
    # Get calculations
    def get_capacitance(radius, length, eh_r):
        ''' assumes a linear wall '''
        return 3.0 * length * np.pi * (radius ** 2) / (2 * eh_r)
    
    def get_inductance(radius, length, density):
        ''' inductance '''
        return length * density / (np.pi * (radius ** 2))
    
    def get_resistance(radius, length, viscosity):
        ''' viscous resistance '''
        return 8.0 * viscosity * length / (np.pi * (radius ** 4))
    
    
    
    mpa = {
            "boundary_conditions": {
                "inlet": "INFLOW"
            },
            "vessel_id": 0,
            "vessel_length": length, # placeholder
            "vessel_name": "mpa",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "C": get_capacitance(avg_radius, length, params.linear_ehr),
                "L": get_inductance(avg_radius, length, params.density),
                "R_poiseuille": get_resistance(avg_radius, length, params.viscosity),
            }
        }
    lpa =  [{
            "vessel_id": 1,
            "vessel_length": 10.0,
            "vessel_name": "lpa0",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "C": 0,
                "R_poiseuille": 0,
            }
        },
            { "boundary_conditions": {
                "outlet": "LPA_BC"
            },
            "vessel_id": 2,
            "vessel_length": 10.0,
            "vessel_name": "lpa1",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }]
    rpa =  [{ 
            "vessel_id": 3,
            "vessel_length": 10.0,
            "vessel_name": "rpa0",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                'C': 0,
                "R_poiseuille": 0,
            }
        }, 
            { "boundary_conditions": {
                "outlet": "RPA_BC"
            },
            "vessel_id": 4,
            "vessel_length": 10.0,
            "vessel_name": "rpa1",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }]
    solver_data.vessel = [mpa] + lpa + rpa
    
    return solver_data
    
def modify_params(solver_data: Solver0D, x ):
    ''' modifies the optimizer variables '''
    vess = solver_data.vessel

    # LPA
    vess[1]['zero_d_element_values']['R_poiseuille']  = x[0]
    vess[1]['zero_d_element_values']['C'] = x[1]
    vess[2]['zero_d_element_values']['R_poiseuille'] = x[2]
    
    # RPA
    vess[3]['zero_d_element_values']['R_poiseuille']  = x[3]
    vess[3]['zero_d_element_values']['C'] =  x[4]
    vess[4]['zero_d_element_values']['R_poiseuille'] = x[5]
    
def get_initial_cond(params: TuneParams, inflow: Inflow, solver_data: Solver0D):
    ''' initial cond with numpy array in form of [Rp_LPA, C_LPA, Rd_LPA, Rp_RPA, C_RPA, Rd_RPA]'''
    PVR = (params.mPAP_meas - params.cap_wedge_pressure)/ inflow.mean_inflow
    mpa_cap = solver_data.vessel[0]['zero_d_element_values']['C']
    x0 = np.array([.4 * PVR,
                    mpa_cap,
                    1.6 * PVR,
                    .4 * PVR,
                    mpa_cap,
                    1.6 * PVR]
                    )
    return x0

################
# Optimization #
################

def termination_closure(params: TuneParams):
        params = {'tolerance': params.tolerance, 'losses': [],'patience': params.pat, 'diff_tol': params.pat_tol, 'cur_patience': 0 }

        def termination_callback(xk, state):
            # terminate if loss has not improved
            if params['losses']:
                if params['losses'][-1] - state['fun'] < params['diff_tol']:
                    if state['cg_stop_cond'] == 1 : # only if iteration limit reached
                        params['cur_patience'] +=1
                else:
                    if state['cg_stop_cond'] != 0: # only if it actually evaluated something
                        params['cur_patience'] = 0
            
            params['losses'].append(state['fun'])
                
            if params['cur_patience'] >= params['patience']:
                print('Optimization exited abnormally: target tolerance not reached.')
                return True
                
            # terminate if loss is below tolerance
            if state['fun'] < params['tolerance']:
                print(f"Optimization exited normally: target tolerance {params['tolerance']} reached.")
                return True
            else: 
                return False
        return termination_callback

def opt_function(x, tune_params: TuneParams, tune_solver, inflow: Inflow):
    
    # run a simulation
    run_sim(x, tune_solver, last_cycle=True)
    
    # extract results/compute loss
    results = np.load(get_result_file(tune_solver), allow_pickle=True).item()
    
    mPAP_loss, Q_RPA_loss, flow_mse, _ = loss_function(results, tune_params, inflow, compute_last_cycle=False)

    return mPAP_loss + Q_RPA_loss + flow_mse
    
def loss_function(results, tune_params: TuneParams, inflow: Inflow, compute_last_cycle = False):
    ''' loss function'''
    
    if compute_last_cycle:
        last_cycle = -1 * tune_params.num_timesteps_per_cycle
        time = results['time'][last_cycle:]
        time = time - time[0] + results['time'][0]
    else:
        last_cycle = 0
        time = results['time']

    mPAP_sim = np.trapz(results['pressure']['P_BC0_inlet_V0'][last_cycle:], time) / inflow.tc
    Q_RPA_sim = np.trapz(results['flow']['Q_V4_BC4_outlet'][last_cycle:], time) / inflow.tc
    mPAP_meas = tune_params.mPAP_meas
    Q_RPA_meas = inflow.mean_inflow * tune_params.rpa_flow_split
    
    f = interp1d(time, results['flow']['Q_V4_BC4_outlet'][last_cycle:])
    
    mse_loss = np.square(np.divide(f(inflow.t[:-1]) - (inflow.Q[:-1] * tune_params.rpa_flow_split), (inflow.Q[:-1] * tune_params.rpa_flow_split))).mean()

    mPAP_loss = ((mPAP_sim - mPAP_meas) / mPAP_meas)**2
    Q_RPA_loss =  ((Q_RPA_sim - Q_RPA_meas) / Q_RPA_meas) ** 2
    return  mPAP_loss , Q_RPA_loss, mse_loss, (mPAP_sim, Q_RPA_sim, mPAP_meas, Q_RPA_meas)

def run_sim(x, solver_file, last_cycle = True ):
    ''' run simulation:
    last_cycle -> True to only save results for last cycle. False for all cycle results'''
    # read solver file
    solver_data = Solver0D()
    solver_data.read_solver_file(solver_file)

    
    # modify conditions
    modify_params(solver_data, x)
    solver_data.write_solver_file(solver_file=solver_file)             
    
    # run simulation
    blockPrint() # silence
    zerod.solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path=solver_file,
                                    last_cycle=last_cycle,
                                    save_results_all=True,
                                    save_results_branch=False,
                                    use_steady_soltns_as_ics = False)
    enablePrint()

def get_result_file(solver_file):
    return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_all_results.npy'


###########
# Results #
###########

def validate_results(tune_params: TuneParams, inflow: Inflow, sim_results, opt_results, results_dir):

    
    
    ## save inflow graph
    fig,ax = plt.subplots(1,1 )
    ax.plot(inflow.t, inflow.Q)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('flow (ml/s)')
    fig.savefig(os.path.join(results_dir , 'inflow.png'))
    
    
    # save some computations as a json
    results_dict = {}
    x = opt_results['x']
    results_dict['Final Params'] = {'Rp_LPA': x[0], 
                                    'C_LPA': x[1],
                                    'Rd_LPA': x[2],
                                    'Rp_RPA': x[3],
                                    'C_RPA': x[4],
                                    'Rd_RPA': x[5]}
    
    results_dict['columns'] = ['Optimized', 'Desired']
    
    mPAP_loss, Q_RPA_loss, flow_mse, (mPAP_sim, Q_RPA_sim, mPAP_meas, Q_RPA_meas) = loss_function(sim_results, tune_params, inflow, compute_last_cycle=True)
    results_dict['mPAP'] = [mPAP_sim, mPAP_meas]
    results_dict['max_pressure'] = [sim_results['pressure']['P_BC0_inlet_V0'].max(), (m2d(18), m2d(25))]
    results_dict['rpa_flow_split'] = [Q_RPA_sim/inflow.mean_inflow, tune_params.rpa_flow_split]
    results_dict['Q_RPA'] = [Q_RPA_sim, Q_RPA_meas]
    
    results_dict['losses'] = {'mPAP_loss': mPAP_loss,
                              'Q_RPA_loss': Q_RPA_loss,
                              'MSE_loss':flow_mse}
    
    
    
    with open(os.path.join(results_dir, 'values.json'), 'w') as json_file:
        json.dump(results_dict, json_file, indent = 4)
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(2, 3, figsize=(30, 20))
    ax = ax.flatten()
    start = -3 * tune_params.num_timesteps_per_cycle
    ax[0].plot(sim_results['time'][start:], sim_results['pressure']['P_BC0_inlet_V0'][start:]/1333.22)
    ax[0].set_title('Inlet Pressure')
    ax[1].plot(sim_results['time'][start:], sim_results['pressure']['P_V2_BC2_outlet'][start:]/1333.22)
    ax[1].set_title('LPA Outlet Pressure')
    ax[2].plot(sim_results['time'][start:], sim_results['pressure']['P_V4_BC4_outlet'][start:]/1333.22)
    ax[2].set_title('RPA Outlet Pressure')
    
    for i in range(3):
        ax[i].set_xlabel('time (s)')
        ax[i].set_ylabel('pressure (mmHg)')
        
    ax[3].plot(sim_results['time'][start:], sim_results['flow']['Q_BC0_inlet_V0'][start:])
    ax[3].set_title('Inlet Flow')
    ax[4].plot(sim_results['time'][start:], sim_results['flow']['Q_V2_BC2_outlet'][start:])
    ax[4].set_title('LPA Outlet Flow')
    ax[5].plot(sim_results['time'][start:], sim_results['flow']['Q_V4_BC4_outlet'][start:])
    ax[5].set_title('RPA Outlet Flow')
    
    for i in range(3, 6):
        ax[i].set_xlabel('time (s)')
        ax[i].set_ylabel('flow (ml/s)')

    fig.savefig(os.path.join(results_dir, 'waveforms.png'))


def sens_test(opt_results, tune_params : TuneParams, inflow: Inflow, solver_file, sensitivity_dir):
    
    x_opt = opt_results['x']
    mapping = {'Rp_LPA': 0, 'C_LPA': 1, 'Rd_LPA': 2, 'Rp_RPA':3, 'C_RPA':4, 'Rd_RPA':5}
    
    for var_name, index in mapping.items():
        
        fig, ax = plt.subplots(3, 1, figsize = (40, 30))
        ax = ax.flatten()
        ax[0].set_xlabel(f'pct change {var_name}')
        ax[1].set_xlabel(f'pct change {var_name}')
        ax[0].set_ylabel(f'mPAP')
        ax[1].set_ylabel(f'Q_RPA_avg')
        ax[0].set_title('mPAP change')
        ax[1].set_title('Q_RPA change')
        ax[2].set_xlabel(f'pct_change {var_name}')
        ax[2].set_ylabel(f'MSE')
        ax[2].set_title('MSE')
        #ax[3].set_xlabel(f'time (s)')
        #ax[3].set_ylabel(f'flow (ml/s)')
        #ax[4].
        
        mod = np.linspace(.9, 1.1, 20)
        mpap = []
        q_rpa = []
        mses = []
        for pct in mod:
            x_test = np.copy(x_opt)
            x_test[index] *= pct
            
            run_sim(x_test, solver_file , last_cycle = True)
        
            # extract results/compute loss
            results = np.load(get_result_file(solver_file), allow_pickle=True).item()
            
            _, _, mse, (mPAP_sim, Q_RPA_sim, _, _) = loss_function(results,tune_params, inflow, compute_last_cycle = False)
            mpap.append(mPAP_sim)# - mPAP_opt)
            q_rpa.append(Q_RPA_sim)# - Q_RPA_opt)
            mses.append(mse)
            

        
        ax[0].plot(mod - 1, mpap)
        ax[1].plot(mod - 1, q_rpa)
        ax[2].plot(mod - 1, mses)
        
        fig.savefig(os.path.join(sensitivity_dir, f'{var_name}_change.png'))


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

def calc_rcrs(results, inlet, area_file, out_dir, Pd):
    
    areas = load_area_file(area_file)
    
    del areas[inlet]
    
    rcrs = split_bc(areas, results['x'])
    
    bcs = BoundaryConditions()
    
    for name, vals in rcrs.items():
        bcs.add_rcr(face_name=name, Rp = vals['Rp'], C = vals['C'], Rd = vals['Rd'], Pd = Pd)
    
    bcs.write_rcrt_file(out_dir)
    
    return rcrs

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

###########
# Drivers #
###########

def tool_main(args):
    raise NotImplementedError()
    # TODO: PARSE A CONFIG FILE for tuning parameters and update Tune Params

def dev_main(args: Namespace):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        print(f'Tuning {model_name}...')
        model = org.find_model(model_name)
        
        # perform optimization
        
        
        # ONLY FOR HEALTHY, perform optimization
        if mod_path.type == 'healthy':
        
            results_file = os.path.join(mod_path.solver_files, model_name + '_results.npy')
            # use defaults
            tuning_params = TuneParams()
            if os.path.exists(results_file):
                print('Result file found. Skipping optimization')
                results = np.load(results_file, allow_pickle = True).item()
            else:
                
                # set up inflow file
                if mod_path.params['files']['rom_inflow']:
                    inflow = Inflow(mod_path.params['files']['rom_inflow'], inverse = False, smooth = True, n_points=1000)
                elif mod_path.params['files']['inflow']:
                    inflow = Inflow(mod_path.params['files']['inflow'], inverse = True, smooth = True, n_points=1000)
                else:
                    raise ValueError('inflow files do not exist')

                # load centerline file
                centerlines = Centerlines()
                centerlines.load_centerlines(mod_path.tune_centerlines)
                
                
                
                # set up solver file
                solver_data = write_0d_dict(tuning_params, inflow, centerlines, mod_path.params['model']['units'])
                x0 = get_initial_cond(tuning_params, inflow, solver_data)
                modify_params(solver_data, x0)
                solver_data.write_solver_file(mod_path.tune_solver)
                
                # run optimizer
                bounds = optimize.Bounds([0,0,0,0,0,0], [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf], keep_feasible=True)
                results = optimize.minimize(opt_function,
                                    x0,
                                    (tuning_params,
                                    mod_path.tune_solver,
                                    inflow),
                                    method='trust-constr',
                                    bounds=bounds,
                                    options={
                                                'verbose':3,
                                                }, # to test
                                    callback=termination_closure(tuning_params))
                
            
                
                np.save(results_file,convert_to_dict(results), allow_pickle=True)
                
                run_sim(results['x'], mod_path.tune_solver, last_cycle=False)
                
                # validate results
                sim_results = np.load(get_result_file(mod_path.tune_solver), allow_pickle = True).item()
                validate_results(tuning_params, inflow, sim_results, results, mod_path.solver_files)
                
                if args.sens_test:
                    print('Running sensitivity tests...', end = '\t')
                    sens_dir = os.path.join(mod_path.solver_files, 'sens_test')
                    if not os.path.exists(sens_dir):
                        os.mkdir(sens_dir)
                    sens_test(results,tuning_params, inflow, mod_path.tune_solver, sens_dir)
                    
                    run_sim(results['x'], mod_path.tune_solver, last_cycle=False)
                    print('Done')
            

            calc_rcrs(results, mod_path.params['model']['inlet'], mod_path.params['files']['cap_info'], mod_path.solver_files, tuning_params.cap_wedge_pressure )
            print('Done')
        
        # Only of NCI_stenosis
        if mod_path.type == 'nci_stenosis':
            # moves already existing rcrt path to appropriate directory and remaps it
            
            if not mod_path.params['files']['rcrt_file']:
                raise ValueError('Missing path to rcrt.dat file')
            
            if args.ingrid:
                
                convert_old_rcrt(inlet=mod_path.params['model']['inlet'],
                                mdl_cvpre_file=mod_path.params['files']['cvpre_file'],
                                old_rcrt_file=mod_path.params['files']['rcrt_file'],
                                out_dir=mod_path.solver_files,
                                solver3d=mod_path.params['files']['solver3d_file'],
                                ingrid=args.ingrid)
            else:
                
                convert_old_rcrt(inlet=mod_path.params['model']['inlet'],
                                mdl_cvpre_file=mod_path.params['files']['mdl_file'],
                                old_rcrt_file=mod_path.params['files']['rcrt_file'],
                                out_dir=mod_path.solver_files,
                                solver3d=mod_path.params['files']['solver3d_file'],
                                ingrid=args.ingrid)
           
                
            
            
        
        
        
        
        
        
        
        
        
        



if __name__ == '__main__':
    parser, dev, tool = create_parser(description = 'Tunes a model and outputs an rcr bc file.')
    
    # dev params
    dev.add_argument('-root', dest = 'root', type = str, default = '.',  help = 'Root to entire project')
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-sens_test', action = 'store_true', default = False, help = 'whether to run sensitivity tests or not')
    #dev.add_argument('-ingrid', action = 'store_true', default = False, help = 'if mapping from one of ingrids models')
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)
    