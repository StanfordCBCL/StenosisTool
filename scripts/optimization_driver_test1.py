


from attr import s
from src.flow import Inflow0D
from src.polydata import Centerlines
from src.file_io import check_exists_bool, parse_face_names
from src.lpn import Solver0D
from src.bc import BoundaryConditions
from src.solver_results import SolverResults
from src.misc import m2d, d2m, create_tool_parser, get_solver_name
from src.run_sim import run_sim

import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
import json
from tqdm import tqdm
from functools import partialmethod
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
        self.mPAP_meas = [m2d(12), m2d(16)] # mmHg -> cgs
        self.cap_wedge_pressure = m2d(7) # mmHg -> cgs
        self.viscosity = .04 # cgs
        self.density = 1.06 # cgs
        self.linear_ehr = 1.2e6 # dyne/cm^2
        self.maxPAP_meas = [m2d(18), m2d(25)]
        
        
        ## for termination
        self.pat = 5
        self.pat_tol  = 1e-6
        self.tolerance = 0.001


############################
# Constructing Tuner Model #
############################
def compute_total_resistance_round1(cur_node: Solver0D.Node):
    ''' computes total resistance coefficients of a tree recursively
    '''
    
    if not cur_node.children:
        r_p = cur_node.vessel_info[0]['zero_d_element_values']['R_poiseuille']
        return r_p
    
    total_inv_res = 0
    for child_node in cur_node.children:
        r_p_child = compute_total_resistance_round1(child_node)
        if r_p_child != 0:
            total_inv_res += 1/r_p_child
    total_child_res = 1/total_inv_res
    r_p = cur_node.vessel_info[0]['zero_d_element_values']['R_poiseuille'] + total_child_res
    return r_p

def compute_lpa_rpa_resistances_round1(dummy_solver: Solver0D ):
    ''' computes LPA and RPA resistances from each branch
    '''
    vessel_tree = dummy_solver.get_vessel_tree()
    assert len(vessel_tree.children) == 2, 'Pulmonary vasculature should only branch into 2 sections, LPA and RPA'
        
    rpa, lpa = get_lpa_rpa(vessel_tree)
    
    lpa_res = compute_total_resistance_round1(lpa)
    rpa_res = compute_total_resistance_round1(rpa)

    return lpa_res, rpa_res

def get_avg_flow(df):
    flow = np.array(df['flow_in'])
    time = np.array(df['time'])
    return np.trapz(flow, time) / (time[-1] - time[0])

def compute_total_resistance_round2(cur_node: Solver0D.Node, solver_results: SolverResults, tc):
    ''' computes total resistance coefficients of a tree recursively INCLUDING stenosis coefficients
    '''
    if not cur_node.children:
        # resistance factor
        r_p = cur_node.vessel_info[0]['zero_d_element_values']['R_poiseuille']
        # R expansion factor
        r_p += cur_node.vessel_info[0]['zero_d_element_values']['stenosis_coefficient'] * abs(get_avg_flow(SolverResults.only_last_cycle(solver_results.vessel_df('V' + str(cur_node.vess_id[0])), tc)))
        return r_p
    
    total_inv_res = 0
    for child_node in cur_node.children:
        r_p_child = compute_total_resistance_round2(child_node, solver_results, tc)
        if r_p_child != 0:
            total_inv_res += 1/r_p_child
    total_child_res = 1/total_inv_res
    
    r_p = cur_node.vessel_info[0]['zero_d_element_values']['R_poiseuille'] + total_child_res
    r_p += cur_node.vessel_info[0]['zero_d_element_values']['stenosis_coefficient'] * abs(get_avg_flow(SolverResults.only_last_cycle(solver_results.vessel_df('V' + str(cur_node.vess_id[0])), tc)))
    return r_p

def compute_lpa_rpa_resistances_round2(dummy_solver: Solver0D, solver_results: SolverResults ):
    ''' computes LPA and RPA resistances from each branch
    '''
    vessel_tree = dummy_solver.get_vessel_tree()
    assert len(vessel_tree.children) == 2, 'Pulmonary vasculature should only branch into 2 sections, LPA and RPA'
        
    rpa, lpa = get_lpa_rpa(vessel_tree)
    
    lpa_res = compute_total_resistance_round2(lpa, solver_results, tc = dummy_solver.inflow.tc)
    rpa_res = compute_total_resistance_round2(rpa, solver_results, tc = dummy_solver.inflow.tc)

    return lpa_res, rpa_res

    
def get_lpa_rpa(vessel_tree: Solver0D.VesselNode):
    # get tree
    
    
    first_branch = vessel_tree.children[0]
    second_branch = vessel_tree.children[1]
    
    # follow the first branch till it finishes to determine which is LPA and which is RPA
    if first_branch.side == 'lpa':
        lpa = first_branch
        rpa = second_branch
    else:
        lpa = second_branch
        rpa = first_branch
    return rpa, lpa


    


def write_0d_dict(params: TuneParams, inflow: Inflow0D, centerlines: Centerlines,dummy_solver: Solver0D) -> Solver0D:
    ''' sets up solver file skeleton'''
    solver_data = Solver0D()
    solver_data.setup_empty_solver()
    
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
            "junction_type": "NORMAL_JUNCTION",
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
    units = dummy_solver.simulation_params['units']
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
    
    lpa_res, rpa_res = compute_lpa_rpa_resistances_round1(dummy_solver=dummy_solver)
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
                "R_poiseuille": rpa_res
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
    
    # check that the results are within reasonable bounds
    tree = dummy_solver.get_vessel_tree()
    tree_res = compute_total_resistance_round1(tree)
    print('delta_P =', d2m(inflow.mean_inflow * tree_res))
    
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
    
def get_initial_cond(params: TuneParams, inflow: Inflow0D, solver_data: Solver0D):
    ''' initial cond with numpy array in form of [Rp_LPA, C_LPA, Rd_LPA, Rp_RPA, C_RPA, Rd_RPA]'''
    PVR = ((params.mPAP_meas[0] + params.mPAP_meas[1])/2 - params.cap_wedge_pressure)/ inflow.mean_inflow
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

def opt_function(x, tune_params: TuneParams, tune_solver: Solver0D, inflow: Inflow0D):
    
    # run a simulation
    rez = run_sim_wrapper(x, tune_solver)
    
    losses, _ = loss_function(rez, tune_params, inflow)

    return sum(losses.values())

def squared_error(target, sim):
    return ((sim - target)/target)**2

def piecewise_error(lower, upper, sim):
    if sim > lower and sim < upper:
        return 0
    if sim < lower:
        return squared_error(lower, sim)
    if sim > upper:
        return squared_error(upper, sim)
        
def one_sided_loss(target, sim, mode = 'lower'):
    if mode == 'lower':
        if sim < target:
            return squared_error(target, sim)
        else:
            return 0
    elif mode == 'upper':
        if sim > target:
            return squared_error(target, sim)
        else:
            return 0
    
def loss_function(results: SolverResults, tune_params: TuneParams, inflow: Inflow0D):
    ''' loss function'''

    mpa = SolverResults.only_last_cycle(results.vessel_df('V0'), inflow.tc)
    rpa = SolverResults.only_last_cycle(results.vessel_df('V6'), inflow.tc)
    time = mpa['time'].to_numpy()
        
    #mPAP
    mPAP_sim = np.trapz(mpa['pressure_in'].to_numpy(), mpa['time'].to_numpy()) / inflow.tc
    mPAP_meas = tune_params.mPAP_meas
    mPAP_loss = piecewise_error(mPAP_meas[0], mPAP_meas[1], mPAP_sim)
    
    #maxPAP
    mpa_pressure = mpa['pressure_in'].to_numpy()
    maxPAP_sim = mpa_pressure.max()
    maxPAP_meas = tune_params.maxPAP_meas
    maxPAP_loss = piecewise_error(maxPAP_meas[0], maxPAP_meas[1], maxPAP_sim)
    
    #maxPAP_t = inflow_max_t
    maxPAP_t = time[np.where(mpa_pressure == maxPAP_sim)][0]
    maxPAP_t_loss = squared_error(inflow.max_inflow_t, maxPAP_t)
    
    #minPAP > PCWP
    minPAP_sim = mpa_pressure.min()
    PCWP = tune_params.cap_wedge_pressure
    minPAP_loss = one_sided_loss(PCWP, minPAP_sim, mode = 'lower')
    
    #minPAP_t = inflow_min_t
    minPAP_t = time[np.where(mpa_pressure == minPAP_sim)][0]
    minPAP_t_loss = squared_error(inflow.min_inflow_t, minPAP_t)
    
    #flow split
    mQ_RPA_sim = np.trapz(rpa['flow_out'].to_numpy(), rpa['time'].to_numpy()) / inflow.tc
    mQ_RPA_meas = inflow.mean_inflow * tune_params.rpa_flow_split
    mQ_RPA_loss =  squared_error(mQ_RPA_meas, mQ_RPA_sim)
    
    # curvature
    f = interp1d(mpa['time'].to_numpy(), rpa['flow_out'].to_numpy())
    mse_loss = np.square(np.divide(f(inflow.t[:-1]) - (inflow.Q[:-1] * tune_params.rpa_flow_split), (inflow.Q[:-1] * tune_params.rpa_flow_split))).mean() 
    
    
    losses = {'mPAP_loss': mPAP_loss,
              'maxPAP_loss': maxPAP_loss,
              'maxPAP_t_loss': maxPAP_t_loss,
              'minPAP_loss': minPAP_loss,
              'minPAP_t_loss': minPAP_t_loss,
              'mQ_RPA_loss': mQ_RPA_loss,
              'mse_loss': 0.01 * mse_loss}

    vals = {'mPAP_sim': mPAP_sim,
            'mPAP_meas': mPAP_meas,
            'maxPAP_sim': maxPAP_sim,
            'maxPAP_meas': maxPAP_meas,
            'maxPAP_t': maxPAP_t,
            'minPAP_sim': minPAP_sim,
            'minPAP_meas': PCWP,
            'minPAP_t': minPAP_t,
            'mQ_RPA_sim': mQ_RPA_sim,
            'mQ_RPA_meas': mQ_RPA_meas}

    return  losses, vals

def run_sim_wrapper(x, solver: Solver0D):
    ''' run simulation:
    last_cycle -> True to only save results for last cycle. False for all cycle results'''

    # modify conditions
    modify_params(solver, x)       
    
    # run simulation
    rez = run_sim(solver=solver,
                  use_steady_soltns=True,
                  save_branch_results=False,
                  save_csv=False,
                  debug=False)
    
    return rez


###########
# Results #
###########

def validate_results(tune_params: TuneParams, inflow: Inflow0D, sim_results: SolverResults, opt_results: optimize.OptimizeResult, results_dir):
    
    
    # save some computations as a json
    results_dict = {}
    x = opt_results['x']
    results_dict['Final Params'] = {'Rp_LPA': x[0], 
                                    'C_LPA': x[1],
                                    'Rd_LPA': x[2],
                                    'Rp_RPA': x[3],
                                    'C_RPA': x[4],
                                    'Rd_RPA': x[5]}
    v0 = sim_results.vessel_df('V0')
    results_dict['columns'] = ['Optimized', 'Desired']
    losses, vals = loss_function(sim_results, tune_params, inflow)
    results_dict['losses'] = losses
    results_dict['sim/targets'] = vals
    
    
    with open(os.path.join(results_dir, 'values.json'), 'w') as json_file:
        json.dump(results_dict, json_file, indent = 4)
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(2, 3, figsize=(30, 20))
    ax = ax.flatten()
    ltc = -3 * tune_params.num_timesteps_per_cycle
    
    
    v3 = sim_results.vessel_df('V3')
    v6 = sim_results.vessel_df('V6')
    
    ax[0].plot(v0['time'].to_numpy()[ltc:], d2m(v0['pressure_in'].to_numpy())[ltc:])
    ax[0].set_title('Inlet Pressure', fontdict={'fontsize': 24})
    ax[1].plot(v3['time'].to_numpy()[ltc:], d2m(v3['pressure_out'].to_numpy())[ltc:])
    ax[1].set_title('LPA Outlet Pressure', fontdict={'fontsize': 24})
    ax[2].plot(v6['time'].to_numpy()[ltc:], d2m(v6['pressure_out'].to_numpy())[ltc:])
    ax[2].set_title('RPA Outlet Pressure', fontdict={'fontsize': 24})
    
    for i in range(3):
        ax[i].tick_params(axis="x", labelsize=16) 
        ax[i].tick_params(axis = 'y', labelsize=16)
        ax[i].set_xlabel('time (s)', fontdict={'fontsize': 20})
        ax[i].set_ylabel('pressure (mmHg)', fontdict={'fontsize': 20})
        
    ax[3].plot(v0['time'].to_numpy()[ltc:], v0['flow_in'].to_numpy()[ltc:])
    ax[3].set_title('Inlet Flow', fontdict={'fontsize': 24})
    ax[4].plot(v3['time'].to_numpy()[ltc:], v3['flow_in'].to_numpy()[ltc:])
    ax[4].set_title('LPA Outlet Flow', fontdict={'fontsize': 24})
    ax[5].plot(v6['time'].to_numpy()[ltc:], v6['flow_in'].to_numpy()[ltc:])
    ax[5].set_title('RPA Outlet Flow', fontdict={'fontsize': 24})
    
    for i in range(3, 6):
        ax[i].tick_params(axis="x", labelsize=16) 
        ax[i].tick_params(axis = 'y', labelsize=16)
        ax[i].set_xlabel('time (s)', fontdict={'fontsize': 20})
        ax[i].set_ylabel('flow (ml/s)', fontdict={'fontsize': 20})

    fig.savefig(os.path.join(results_dir, 'waveforms.png'))

####################
# Sensitivity Test #
####################
'''
def sens_test(x, tune_params : TuneParams, inflow: Inflow0D, solver: Solver0D, sensitivity_dir):
    mapping = {'Rp_LPA': 0, 'C_LPA': 1, 'Rd_LPA': 2, 'Rp_RPA':3, 'C_RPA':4, 'Rd_RPA':5}
    
    for var_name, index in mapping.items():
        print(f'\tRunning {var_name} sensitivity test...', end = '\t', flush = True)
        
        fig, ax = plt.subplots(2, 2, figsize = (40, 30))
        ax = ax.flatten()
        ax[0].set_xlabel(f'pct change {var_name}', fontdict={'fontsize': 20})
        ax[1].set_xlabel(f'pct change {var_name}', fontdict={'fontsize': 20})
        ax[0].set_ylabel(f'mPAP', fontdict={'fontsize': 20})
        ax[1].set_ylabel(f'Q_RPA_avg', fontdict={'fontsize': 20})
        ax[0].set_title('mPAP change', fontdict={'fontsize': 24})
        ax[1].set_title('Q_RPA change', fontdict={'fontsize': 24})
        ax[2].set_xlabel(f'pct_change {var_name}', fontdict={'fontsize': 20})
        ax[2].set_ylabel(f'MSE', fontdict={'fontsize': 20})
        ax[2].set_title('MSE', fontdict={'fontsize': 24})
        ax[3].set_xlabel(f'pct_change {var_name}', fontdict = {'fontsize':20})
        ax[3].set_ylabel(f'maxPAP', fontdict={'fontsize': 20})
        ax[3].set_title('maxPAP', fontdict={'fontsize': 24})
        
        for i in range(3):
            ax[i].tick_params(axis="x", labelsize=16) 
            ax[i].tick_params(axis = 'y', labelsize=16)
        
        mod = np.linspace(.8, 1.2, 40)
        mpap = []
        q_rpa = []
        mses = []
        maxpap = []
        for pct in mod:
            x_test = np.copy(x)
            x_test[index] *= pct
            
            rez = run_sim_wrapper(x_test, solver)
            
            losses = (mPAP_sim, Q_RPA_sim, maxPAP_sim, _, _,_) = loss_function(rez, tune_params, inflow)
            mpap.append(mPAP_sim)
            q_rpa.append(Q_RPA_sim)
            mses.append(mse)
            maxpap.append(maxPAP_sim)
            

        
        ax[0].plot(mod - 1, mpap)
        ax[1].plot(mod - 1, q_rpa)
        ax[2].plot(mod - 1, mses)
        ax[3].plot(mod - 1, maxpap)
        
        fig.savefig(os.path.join(sensitivity_dir, f'{var_name}_change.png'))
        print('Done')
'''
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
    
    return bcs

def add_rcrs(solver: Solver0D, bcs: BoundaryConditions):
    ''' Adds appropriate rcrs to the dummy solver 
    '''
    
    bc_inv_map = bcs.get_bc_map()
    
    for bc in solver.bc:
        if bc['bc_type'] == 'RCR':
            bc_values = bc_inv_map[solver.bc_map[bc['bc_name']]]
            bc['bc_values']['Rp'] = bc_values['Rp']
            bc['bc_values']['Rd'] = bc_values['Rd']
            bc['bc_values']['C'] = bc_values['C']
            bc['bc_values']['Pd'] = bc_values['Pd']
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

def main(args: Namespace):
    
    for solver_dir in args.solver_dirs:
         # load full model dummy solver
        join = lambda file: os.path.join(solver_dir, file)
        
        dummy_solver = Solver0D()
        dummy_solver.read_solver_file(join(get_solver_name(solver_dir)))
        print(f"Optimizing model {dummy_solver.simulation_params['model_name']}...")
        
        
        tuning_dir = join('tuning_dir')
        if not os.path.exists(tuning_dir):
            os.mkdir(tuning_dir)
            
        join_tune = lambda file: os.path.join(tuning_dir, file)
        results_file = join_tune('optimize_results.npy')
        if not args.force and check_exists_bool(results_file,ignore = True):
                print(f'Result file found. Skipping optimization for {solver_dir}')
                continue
        
        #! If config is to be implemented
        # if args.config:
        #
        tuning_params = TuneParams()
        
        inflow_file = join('inflow.flow')
        inflow = Inflow0D.from_file(inflow_file)

        # load centerlines file
        centerlines = Centerlines()
        centerlines.load_centerlines(join('model_centerlines.vtp'))
                
        # set up solver file
        tuning_data = write_0d_dict(tuning_params, inflow, centerlines, dummy_solver)
        x0 = get_initial_cond(tuning_params, inflow, tuning_data)
        modify_params(tuning_data, x0)
        
        # run optimizer round 1
        bounds = optimize.Bounds([0,0,0,0,0,0], [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf], keep_feasible=True)
        results = optimize.minimize(opt_function,
                            x0,
                            (tuning_params,
                            tuning_data,
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
        
        # save int
        if args.int:
            print('Saving intermediate results...', end = '\t', flush = True)
            sim_results = run_sim_wrapper(results['x'], tuning_data)
            int_dir = join_tune('int_dir')
            if not os.path.exists(int_dir):
                os.mkdir(int_dir)
            int_solver = os.path.join(int_dir, 'tune_solver.in')
            tuning_data.write_solver_file(int_solver)
            validate_results(tuning_params, inflow, sim_results, results, int_dir)
            print('Done')
            
        # split rcrs
        area_file = join('CapInfo')
        inlet = parse_face_names(join('inlet_face_names.dat'))[0]
        bcs = calc_rcrs(results=results_dict, 
                  inlet=inlet,
                  area_file=area_file,
                  out_dir=solver_dir)
        
        # add rcrs to dummy solver
        add_rcrs(dummy_solver, bcs)
        dummy_solver.write_solver_file(join(get_solver_name(solver_dir)))
        
        # run a simulation
        dummy_rez = run_sim(dummy_solver, use_steady_soltns = True, debug = True)
        
        # compute new tuning model resistance
        lpa_r, rpa_r = compute_lpa_rpa_resistances_round2(dummy_solver, dummy_rez)
        print('New dP:', d2m(inflow.mean_inflow * 1/(1/lpa_r + 1/rpa_r)), 'LPA_R:',lpa_r, 'RPA_R:',rpa_r)
        tuning_data._construct_vessel_map()
        lpa = tuning_data.get_vessel(1)
        lpa['zero_d_element_values']['R_poiseuille'] =  lpa_r
        rpa = tuning_data.get_vessel(4)
        rpa['zero_d_element_values']['R_poiseuille'] =  rpa_r
        
        # use previous x as initial condition
        x0 = results['x']
        modify_params(tuning_data, x0)
        
        # change flow split tuning params
        '''
        vessel_tree = dummy_solver.get_vessel_tree()
        rpa, lpa = get_lpa_rpa(vessel_tree)
        lpa_df = SolverResults.only_last_cycle(dummy_rez.vessel_df('V' + str(lpa.vess_id[-1])), dummy_solver.inflow.tc)
        rpa_df = SolverResults.only_last_cycle(dummy_rez.vessel_df('V' + str(rpa.vess_id[-1])), dummy_solver.inflow.tc)
        time = np.array(lpa_df['time'])
        lpa_flow = np.trapz(np.array(lpa_df['flow_in']), time)/dummy_solver.inflow.tc
        rpa_flow = np.trapz(np.array(rpa_df['flow_in']), time)/dummy_solver.inflow.tc
        rpa_flow_split = rpa_flow/(rpa_flow + lpa_flow)
        print('New RPA Split:', rpa_flow_split)
        tuning_params.rpa_flow_split = rpa_flow_split
        '''
        
        # rerun optimizer
        bounds = optimize.Bounds([0,0,0,0,0,0], [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf], keep_feasible=True)
        results = optimize.minimize(opt_function,
                            x0,
                            (tuning_params,
                            tuning_data,
                            inflow),
                            method='trust-constr',
                            bounds=bounds,
                            options={
                                        'verbose':3,
                                        }, # to test
                            callback=termination_closure(tuning_params))
        
        # validate results
        sim_results = run_sim_wrapper(results['x'], tuning_data)
        tune_solver = join_tune('tune_solver.in')
        tuning_data.write_solver_file(tune_solver)
        validate_results(tuning_params, inflow, sim_results, results, tuning_dir)

        '''
        # run sensitivity test
        if args.sens_test:
            print('Running sensitivity tests...')
            sens_dir = join_tune('sens_test')
            if not os.path.exists(sens_dir):
                os.mkdir(sens_dir)
            sens_test(results['x'], tuning_params, inflow, tuning_data, sens_dir)
            print('Done')
        '''
        
        # repeat
        # split rcrs
        area_file = join('CapInfo')
        inlet = parse_face_names(join('inlet_face_names.dat'))[0]
        bcs = calc_rcrs(results=results_dict, 
                  inlet=inlet,
                  area_file=area_file,
                  out_dir=solver_dir)
        
        # add rcrs to dummy solver
        add_rcrs(dummy_solver, bcs)
        dummy_solver.write_solver_file(join(get_solver_name(solver_dir)))
        
        print('Done')

        
        



if __name__ == '__main__':

    tool = create_tool_parser(desc = 'Tunes a model')
    
    # dev params
    tool.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Whether to run an optimization if one already exists')
    tool.add_argument('-int', action = 'store_true', default = False, help = 'intermediate results')
    tool.add_argument('-sens_test', action = 'store_true', default = False, help = 'whether to run sensitivity tests or not')
    args = tool.parse_args()
    
    
    main(args)
    