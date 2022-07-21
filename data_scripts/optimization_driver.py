from src.data_org import DataPath
from src.flow import Inflow
from src.centerlines import Centerlines
from src.file_io import Solver0D, check_exists_bool
from src.run_sim import get_result_file
from src.misc import m2d, d2m, create_parser
from src.run_sim import run_sim, get_result_file

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

##########
# Params #
##########

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


############################
# Constructing Tuner Model #
############################
def compute_total_resistance(cur_node: Solver0D.Node):
    ''' computes total resistance of a tree recursively'''
    if not cur_node.children:
        return cur_node.vessel_info['zero_d_element_values']['R_poiseuille']
    
    total_inv_res = 0
    for child_node in cur_node.children:
        total_inv_res += 1/compute_total_resistance(child_node)
    total_child_res = 1/total_inv_res
    
    return cur_node.vessel_info['zero_d_element_values']['R_poiseuille'] + total_child_res
    
        

def compute_lpa_rpa_resistances(dummy_solver: Solver0D ):
    ''' computes LPA and RPA resistances from each branch'''
    
    # map rcrt bc to a name
    outlet_face_names_file = os.path.join(os.path.dirname(dummy_solver.solver_file), 'outlet_face_names.dat')
    with open(outlet_face_names_file, 'r') as ofile:
        outlet_face_names = ofile.readlines()
    rcrt_map = {}
    for i, name in enumerate(outlet_face_names):
        rcrt_map['RCR_' + str(i)] = name.rstrip()
        
    # get tree
    vessel_tree = dummy_solver.get_vessel_tree()
    assert len(vessel_tree.children) == 2, 'Pulmonary vasculature should only branch into 2 sections, LPA and RPA'
    
    first_branch = vessel_tree.children[0]
    second_branch = vessel_tree.children[1]
    
    # follow the first branch till it finishes to determine which is LPA and which is RPA
    cur_branch = first_branch
    while cur_branch.children:
        cur_branch = cur_branch.children[0]
    rcr_name = rcrt_map[cur_branch.vessel_info['boundary_conditions']['outlet']]
    if 'lpa' in rcr_name.lower():
        lpa = first_branch
        rpa = second_branch
    else:
        lpa = second_branch
        rpa = first_branch
    
    lpa_res = compute_total_resistance(lpa)
    rpa_res = compute_total_resistance(rpa)
    
    return lpa_res, rpa_res
    
    
    


def write_0d_dict(params: TuneParams, inflow: Inflow, centerlines: Centerlines,dummy_solver: Solver0D, units = 'mm') -> Solver0D:
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
                4
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
                2
            ],
            "junction_name": "J2",
            "junction_type": "internal_junction",
            "outlet_vessels": [
                3
            ]
        },
        {
            "inlet_vessels": [
                4
            ],
            "junction_name": "J3",
            "junction_type": "internal_junction",
            "outlet_vessels": [
                5
            ]
        },
        {
            "inlet_vessels": [
                5
            ],
            "junction_name": "J4",
            "junction_type": "internal_junction",
            "outlet_vessels": [
                6
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
    
    lpa_res, rpa_res = compute_lpa_rpa_resistances(dummy_solver=dummy_solver)
    lpa =  [{
            "vessel_id": 1,
            "vessel_length": 10.0,
            "vessel_name": "lpa0",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": lpa_res,
            }
        },
            {
            "vessel_id": 2,
            "vessel_length": 10.0,
            "vessel_name": "lpa1",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "C": 0,
                "R_poiseuille": 0,
            }
        },
            { "boundary_conditions": {
                "outlet": "LPA_BC"
            },
            "vessel_id": 3,
            "vessel_length": 10.0,
            "vessel_name": "lpa2",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }]
    rpa =  [{ 
            "vessel_id": 4,
            "vessel_length": 10.0,
            "vessel_name": "rpa0",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": rpa_res,
            }
        }, 
            { 
            "vessel_id": 5,
            "vessel_length": 10.0,
            "vessel_name": "rpa1",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                'C': 0,
                "R_poiseuille": 0,
            }
        }, 
            { "boundary_conditions": {
                "outlet": "RPA_BC"
            },
            "vessel_id": 6,
            "vessel_length": 10.0,
            "vessel_name": "rpa2",
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
    vess[2]['zero_d_element_values']['R_poiseuille']  = x[0]
    vess[2]['zero_d_element_values']['C'] = x[1]
    vess[3]['zero_d_element_values']['R_poiseuille'] = x[2]
    
    # RPA
    vess[5]['zero_d_element_values']['R_poiseuille']  = x[3]
    vess[5]['zero_d_element_values']['C'] = x[4]
    vess[6]['zero_d_element_values']['R_poiseuille'] = x[5]
    
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
    run_sim_wrapper(x, tune_solver, last_cycle=True)
    
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
    Q_RPA_sim = np.trapz(results['flow']['Q_V6_BC6_outlet'][last_cycle:], time) / inflow.tc
    mPAP_meas = tune_params.mPAP_meas
    Q_RPA_meas = inflow.mean_inflow * tune_params.rpa_flow_split
    
    f = interp1d(time, results['flow']['Q_V6_BC6_outlet'][last_cycle:])
    
    mse_loss = np.square(np.divide(f(inflow.t[:-1]) - (inflow.Q[:-1] * tune_params.rpa_flow_split), (inflow.Q[:-1] * tune_params.rpa_flow_split))).mean()

    mPAP_loss = ((mPAP_sim - mPAP_meas) / mPAP_meas)**2
    Q_RPA_loss =  ((Q_RPA_sim - Q_RPA_meas) / Q_RPA_meas) ** 2
    return  mPAP_loss , Q_RPA_loss, mse_loss, (mPAP_sim, Q_RPA_sim, mPAP_meas, Q_RPA_meas)

def run_sim_wrapper(x, solver_file, last_cycle = True ):
    ''' run simulation:
    last_cycle -> True to only save results for last cycle. False for all cycle results'''
    # read solver file
    solver_data = Solver0D()
    solver_data.read_solver_file(solver_file)

    
    # modify conditions
    modify_params(solver_data, x)
    solver_data.write_solver_file(solver_file=solver_file)             
    
    # run simulation
    run_sim(solver_file, save_branches = False, block_print=True, use_steady_soltns=False, last_cycle = last_cycle)


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
    ax[1].plot(sim_results['time'][start:], sim_results['pressure']['P_V3_BC3_outlet'][start:]/1333.22)
    ax[1].set_title('LPA Outlet Pressure')
    ax[2].plot(sim_results['time'][start:], sim_results['pressure']['P_V6_BC6_outlet'][start:]/1333.22)
    ax[2].set_title('RPA Outlet Pressure')
    
    for i in range(3):
        ax[i].set_xlabel('time (s)')
        ax[i].set_ylabel('pressure (mmHg)')
        
    ax[3].plot(sim_results['time'][start:], sim_results['flow']['Q_BC0_inlet_V0'][start:])
    ax[3].set_title('Inlet Flow')
    ax[4].plot(sim_results['time'][start:], sim_results['flow']['Q_V3_BC3_outlet'][start:])
    ax[4].set_title('LPA Outlet Flow')
    ax[5].plot(sim_results['time'][start:], sim_results['flow']['Q_V6_BC6_outlet'][start:])
    ax[5].set_title('RPA Outlet Flow')
    
    for i in range(3, 6):
        ax[i].set_xlabel('time (s)')
        ax[i].set_ylabel('flow (ml/s)')

    fig.savefig(os.path.join(results_dir, 'waveforms.png'))

####################
# Sensitivity Test #
####################
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
            
            run_sim_wrapper(x_test, solver_file , last_cycle = True)
        
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

########
# Misc #
########

def convert_to_dict(opt_results: optimize.OptimizeResult):
    rez = {}
    for key, val in opt_results.items():
        rez[key] = val
    return rez

###########
# Drivers #
###########

def tool_main(args):
    raise NotImplementedError()
    # TODO: PARSE A CONFIG FILE for tuning parameters and update Tune Params
    # TODO: CREATE A MODEL FILE

def dev_main(args: Namespace):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        print(f'Tuning {model_name}...')
        model = org.find_model(model_name)
        
        # perform optimization
        results_file = os.path.join(os.path.dirname(model.tune_solver), model.info['metadata']['name'] + '_optimize_results.npy')
        if not args.force and check_exists_bool(results_file,ignore = True):
            print(f'Result file found. Skipping optimization for {model_name}')
            continue
        
        # use defaults for tuning if config file is not found
        if 'tune_config' in model.info['files']:
            print('Reading From a config file has not yet been implemented. Running defaults.')
            #! Change this when implementing config
            tuning_params = TuneParams()
        else:
            tuning_params = TuneParams()
        
        # set up inflow file
        if model.info['files']['rom_inflow']:
            inflow = Inflow(model.info['files']['rom_inflow'], inverse = False, smooth = True, n_points=1000)
        elif model.info['files']['inflow']:
            inflow = Inflow(model.info['files']['inflow'], inverse = True, smooth = True, n_points=1000)
        else:
            raise ValueError('Inflow files do not exist')

        # load centerlines file
        centerlines = Centerlines()
        centerlines.load_centerlines(model.model_centerlines)
        
        # load full model dummy solver
        dummy_solver = Solver0D()
        dummy_solver.read_solver_file(model.model_solver)
                
        # set up solver file
        solver_data = write_0d_dict(tuning_params, inflow, centerlines, dummy_solver, model.info['model']['units'])
        x0 = get_initial_cond(tuning_params, inflow, solver_data)
        modify_params(solver_data, x0)
        solver_data.write_solver_file(model.tune_solver)
        # run optimizer
        bounds = optimize.Bounds([0,0,0,0,0,0], [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf], keep_feasible=True)
        results = optimize.minimize(opt_function,
                            x0,
                            (tuning_params,
                            model.tune_solver,
                            inflow),
                            method='trust-constr',
                            bounds=bounds,
                            options={
                                        'verbose':3,
                                        }, # to test
                            callback=termination_closure(tuning_params))
        
        # add Pd to results 
        results_dict = convert_to_dict(results)
        results_dict['Pd'] = tuning_params.cap_wedge_pressure
        np.save(results_file, results_dict , allow_pickle=True)
        
        run_sim_wrapper(results['x'], model.tune_solver, last_cycle=False)
        
        # validate results
        sim_results = np.load(get_result_file(model.tune_solver), allow_pickle = True).item()
        validate_results(tuning_params, inflow, sim_results, results, model.tuning_dir)
        
        if args.sens_test:
            print('Running sensitivity tests...', end = '\t')
            sens_dir = os.path.join(model.tuning_dir, 'sens_test')
            if not os.path.exists(sens_dir):
                os.mkdir(sens_dir)
            sens_test(results, tuning_params, inflow, model.tune_solver, sens_dir)
            
            run_sim_wrapper(results['x'], model.tune_solver, last_cycle=False)
        print('Done')

        
        
        
        
        
        

        
        
        
        



if __name__ == '__main__':

    parser, dev, tool = create_parser(desc = 'Tunes a model')
    
    # dev params
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Whether to run an optimization if one already exists')
    dev.add_argument('-c', dest = 'c', action = 'store_true', default = False, help = 'Whether to use c solver')
    dev.add_argument('-sens_test', action = 'store_true', default = False, help = 'whether to run sensitivity tests or not')
    #dev.add_argument('-ingrid', action = 'store_true', default = False, help = 'if mapping from one of ingrids models')
    args = parser.parse_args()
    
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        if args.c:
            print('C solver has not been implemented yet.')
            #! IMPLEMENT C
        dev_main(args)
    