# File: tune_bc.py
# File Created: Monday, 31st October 2022 8:46:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:16:31 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Tunes Boundary Conditions for a 0D model using a simplified nonlinear tuning model.


from svinterface.core.zerod.solver import SolverResults, Solver0Dcpp
from svinterface.core.zerod.lpn import LPN
from svinterface.core.bc import Inflow, RCR
from svinterface.manager import Manager
from svinterface.utils.misc import m2d, d2m
from svinterface.utils.io import write_json

from pathlib import Path
import numpy as np
import argparse
import shutil
from copy import deepcopy
from scipy import optimize
import matplotlib.pyplot as plt


##########
# Params #
##########


class TuneParams():
    """ Tuning Parameter Defaults. """
    
    def __init__(self):
        ## optimizer targets defaults for a healthy model
        self.rpa_flow_split = .55
        self.mPAP_meas = [m2d(12), m2d(16)] # mmHg -> cgs
        self.maxPAP_meas = [m2d(18), m2d(25)]
        self.minPAP_meas = [m2d(8), m2d(12)]

        ## assume wedge pressure
        self.cap_wedge_pressure = m2d(7) # mmHg -> cgs
        
        ## for tracking
        self.iter = 0
        self.full_sim_res = None
        
        ## for RPA & LPA resistances
        self.r_RPA = 0
        self.r_LPA = 0

##############
# Tuning LPN #
##############

def construct_tuning_lpn(params: TuneParams, main_lpn: LPN):
    ''' sets up tuning LPN skeleton in C++ form
    '''
    tuning_lpn = LPN()
    tuning_lpn.setup_empty_lpn()
    
    # load in inflow from main_lpn
    inflow = main_lpn.inflow
    
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
    
    tuning_lpn.bc = [inflow_bc, lpa_bc, rpa_bc]
    
    ## simulation_parameters
    tuning_lpn.simulation_params = deepcopy(main_lpn.simulation_params)
    
    ## junction basic
    tuning_lpn.junctions= [
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
    
    ## dummy mpa_junctions
    mpa = [{
            "boundary_conditions": {
                "inlet": "INFLOW"
            },
            "vessel_id": 0,
            "vessel_length": 10.0,
            "vessel_name": "branch_mpa",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }
           ]
    

    lpa =  [{
            "vessel_id": 1,
            "vessel_length": 10.0,
            "vessel_name": "branch_lpa_tree",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
                "stenosis_coefficient": params.r_LPA,
            }
        },
            {
            "vessel_id": 2,
            "vessel_length": 10.0,
            "vessel_name": "branch_lpa_rp_c",
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
            "vessel_name": "branch_lpa_rd",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }]
    rpa =  [{ 
            "vessel_id": 4,
            "vessel_length": 10.0,
            "vessel_name": "branch_rpa_tree",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
                "stenosis_coefficient": params.r_RPA
                
            }
        }, 
            { 
            "vessel_id": 5,
            "vessel_length": 10.0,
            "vessel_name": "branch_rpa_rp_c",
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
            "vessel_name": "branch_rpa_rd",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 0,
            }
        }]
    tuning_lpn.vessel = mpa + lpa + rpa
    
    return tuning_lpn

############
# Opt Func #
############

def opt_function(x, main_lpn: LPN, tuning_lpn: LPN, params: TuneParams):
    ''' Each iteration of optimization runs this
    '''
    
    # run a simulation
    modify_params(tuning_lpn, x)
    solver = Solver0Dcpp(tuning_lpn,
                         use_steady=False,
                         last_cycle_only=True)
    rez = solver.run_sim()
    
    # compute decomposed loss
    mPAP_loss, qRPA_loss, maxPAP_loss, minPAP_loss = loss_function(results = rez,
                         tune_params = params,
                         inflow = main_lpn.inflow,
                         )
    # aggregate
    loss = mPAP_loss + qRPA_loss + maxPAP_loss + minPAP_loss
    
    # report
    params.iter += 1
    print(f"{params.iter:^15}|{mPAP_loss:^15.5f}|{maxPAP_loss:^15.5f}|{minPAP_loss:^15.5f}|{qRPA_loss:^15.5f}|{loss:^15f}")
    
    return loss

def squared_error(target, sim):
    ''' returns squared error
    '''
    return ((sim - target)/target)**2

def piecewise_error(lower, upper, sim):
    ''' returns a piecewise error
    '''
    if sim > lower and sim < upper:
        return 0
    if sim < lower:
        return squared_error(lower, sim)
    if sim > upper:
        return squared_error(upper, sim)
        
def loss_function(results: SolverResults, tune_params: TuneParams, inflow: Inflow, intermediate = False):
    ''' Loss function

        Tunes mPAP, RPA Flow, systolic PAP, .01 of Q mse in rpa.
        returns mean of 4 losses.
    '''
    
    mpa = results.vessel_df('branch_mpa')
    rpa = results.vessel_df('branch_rpa_rd')
    
    
    # qRPA
    qRPA_sim = np.trapz(rpa['flow_out'].to_numpy(), rpa['time'].to_numpy()) / inflow.tc
    qRPA_meas = inflow.mean_inflow * tune_params.rpa_flow_split
    qRPA_loss = squared_error(qRPA_meas, qRPA_sim )

    # maxPAP
    maxPAP_sim = mpa['pressure_in'].to_numpy().max()
    maxPAP_loss = piecewise_error(tune_params.maxPAP_meas[0], tune_params.maxPAP_meas[1], maxPAP_sim)
    
    # minPAP
    minPAP_sim = mpa['pressure_in'].to_numpy().min()
    minPAP_loss = piecewise_error(tune_params.minPAP_meas[0], tune_params.minPAP_meas[1], minPAP_sim)
    
    # mPAP
    mPAP_sim = np.trapz(mpa['pressure_in'].to_numpy(), mpa['time'].to_numpy()) / inflow.tc
    mPAP_loss = piecewise_error(tune_params.mPAP_meas[0], tune_params.mPAP_meas[1], mPAP_sim )
    #mPAP_loss = squared_error(1/3 * maxPAP_sim + 2/3 * minPAP_sim, mPAP_sim)

    if intermediate:
        return mPAP_loss , qRPA_loss, maxPAP_loss, minPAP_loss, (mPAP_sim, qRPA_sim, maxPAP_sim, minPAP_sim)
    
    return mPAP_loss, qRPA_loss, maxPAP_loss, minPAP_loss

#############
# Split RCR #
#############

def split_rcrs(TM: Manager, x, PCWP):
    ''' splits rcr driver '''
    capinfo = load_area_file(str(TM['workspace']['capinfo']))
    del capinfo[TM['metadata']['inlet']]
    
    rcrs = split_bc(capinfo, x)
    
    bcs = RCR()
    
    for name, vals in rcrs.items():
        bcs.add_rcr(face_name=name, Rp = vals['Rp'], C = vals['C'], Rd = vals['Rd'], Pd = PCWP)
    return bcs
    
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

def validate_caps(areas):
    ''' confirm that all caps have either lpa or rpa in them'''
    names = list(areas.keys())
    for name in names:
        if 'lpa' not in name.lower() and 'rpa' not in name.lower():
            raise ValueError('Unable to identify RPA vs. LPA caps: please rename caps to contain lpa or rpa')
    return 

def load_area_file(area_filepath):
    ''' loads a capinfo file
    '''
    with open(area_filepath, 'r') as afile:
        areas = {}
        afile.readline() # ignore first comment line
        for line in afile:
            line = line.rstrip().split()
            areas[line[0]] = float(line[1])
    return areas

def total_area(areas):
    ''' total area
    '''
    return sum(list(areas.values()))

    
###############
# Tune Driver #
###############

def tune(TM: Manager, main_lpn: LPN, tuning_lpn: LPN, params: TuneParams, tuning_dir: Path):
    ''' Tuning Steps
    '''
    # setup initial conditions
    x0 = get_initial_cond(params, main_lpn, tuning_lpn)
    modify_params(tuning_lpn, x0)
    tuning_lpn_file = tuning_dir / (TM['metadata']['model_name'] + '_tuning.in')
    tuning_lpn.write_lpn_file(str(tuning_lpn_file))
    
    # bounds
    bounds = optimize.Bounds([0,0,0], [ np.inf, np.inf, np.inf], keep_feasible=True)
    
    # run optimizer
    print(f"{'Iteration':^15}|{'mPAP Loss':^15}|{'maxPAP Loss':^15}|{'minPAP Loss':^15}|{'qRPA Loss':^15}|{'Total Loss':^15}")
    print("-" * 15 * 7 + "-" * 6)
    
    results = optimize.minimize(fun = opt_function,
                                x0 = x0,
                                args = (main_lpn,
                                        tuning_lpn, params),
                                method='Nelder-Mead',
                                bounds=bounds,
                                options = {'disp': True})
            
    # set prev x0 as start point
    x0 = results.x
    # split using 1:9 ratio
    x_new = np.array([.1*x0[1],x0[0], .9*x0[1],.1*x0[2],x0[0], .9*x0[2] ])
        
    # split rcrs
    bcs = split_rcrs(TM, x_new, params.cap_wedge_pressure)
    bcs.write_rcrt_file(TM['workspace']['lpn_dir'])
    
    main_lpn.update_rcrs(bcs)
    # update on disk
    main_lpn.update()
    

    # store final
    np.save(str(tuning_dir / 'results.npy'), results.x, allow_pickle=True)
    # save tuning file
    modify_params(tuning_lpn, results.x)
    tuning_lpn.write_lpn_file(str(tuning_lpn_file))
    
    # validate results
    validate_results(params,
                        main_lpn,
                        tuning_lpn,
                        x_new,
                        tuning_dir,
                        TM['tune_params'])
    
    print(x_new)
    print("Optimization Complete.")
    
    

def get_initial_cond(params: TuneParams, main_lpn: LPN, tuning_lpn: LPN):
    ''' initial cond with numpy array in form of [Rp_LPA, C_LPA, Rd_LPA, Rp_RPA, C_RPA, Rd_RPA]'''
    # determine PVR
    PVR = ((params.mPAP_meas[0] + params.mPAP_meas[1])/2 - params.cap_wedge_pressure)/ main_lpn.inflow.mean_inflow
    
    x0 = np.array([
                    1e-3,
                    PVR,
                    PVR]
                    )
    return x0

def modify_params(lpn: LPN, x ):
    ''' modifies the optimizer variables '''
    vess = lpn.vessel
    
    # LPA
    vess[2]['zero_d_element_values']['R_poiseuille'] = .1 * x[1]
    vess[2]['zero_d_element_values']['C'] = x[0]
    vess[3]['zero_d_element_values']['R_poiseuille'] = .9 * x[1]
    
    # RPA
    vess[5]['zero_d_element_values']['R_poiseuille'] = .1 * x[2]
    vess[5]['zero_d_element_values']['C'] = x[0]
    vess[6]['zero_d_element_values']['R_poiseuille'] = .9 * x[2]

def convert_to_dict(opt_results: optimize.OptimizeResult):
    """Converts optimizer results to a dict
    """
    rez = {}
    for key, val in opt_results.items():
        rez[key] = val
    return rez

###########
# Results #
###########

def validate_results(tune_params: TuneParams, main_lpn: LPN, tuning_lpn: LPN, x, save_dir: Path, targets):
    
    
    # save some computations as a json
    results_dict = {}
    results_dict['x'] = {'Rp_LPA': x[0], 
                                    'C_LPA': x[1],
                                    'Rd_LPA': x[2],
                                    'Rp_RPA': x[3],
                                    'C_RPA': x[4],
                                    'Rd_RPA': x[5]}
    
    results_dict['columns'] = ['Optimized', 'Desired']
    
    # run a sim w/ last cycle only
    solver = Solver0Dcpp(tuning_lpn,
                         use_steady=False,
                         last_cycle_only=True)
    tune_results = solver.run_sim()
    
    mPAP_loss, qRPA_loss, maxPAP_loss, minPAP_loss, (mPAP_sim, qRPA_sim, maxPAP_sim, minPAP_sim) = loss_function(tune_results, tune_params, main_lpn.inflow, intermediate = True)
    
    results_dict['mPAP'] = [d2m(mPAP_sim), [d2m(tune_params.mPAP_meas[0]), d2m(tune_params.mPAP_meas[1])]] #d2m(2/3 * minPAP_sim + 1/3 * maxPAP_sim)]
    
    results_dict['max_pressure'] = [d2m(maxPAP_sim), [d2m(tune_params.maxPAP_meas[0]), d2m(tune_params.maxPAP_meas[1])]]
    results_dict['min_pressure'] = [d2m(minPAP_sim), [d2m(tune_params.minPAP_meas[0]), d2m(tune_params.minPAP_meas[1])]]
    
    results_dict['rpa_flow_split'] = [qRPA_sim/main_lpn.inflow.mean_inflow, tune_params.rpa_flow_split]
    
    results_dict['qRPA'] = [qRPA_sim, tune_params.rpa_flow_split * main_lpn.inflow.mean_inflow]
    
    
    
    results_dict['losses'] = {'mPAP_loss': mPAP_loss,
                              'qRPA_loss': qRPA_loss,
                              'minPAP_loss':minPAP_loss,
                              'maxPAP_loss': maxPAP_loss}
    
    write_json(save_dir / 'values.json', results_dict, sort_keys = False)
    
    # run a sim without last cycle
    solver = Solver0Dcpp(tuning_lpn,
                         use_steady=False,
                         last_cycle_only=False)
    tune_results = solver.run_sim()
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(2, 3, figsize=(30, 20))
    ax = ax.flatten()
    ltc = -3 * main_lpn.simulation_params['number_of_time_pts_per_cardiac_cycle']
    
    mpa = tune_results.vessel_df('branch_mpa')
    lpa = tune_results.vessel_df('branch_lpa_rd')
    rpa = tune_results.vessel_df('branch_rpa_rd')
    
    ax[0].plot(mpa['time'].to_numpy()[ltc:], d2m(mpa['pressure_in'].to_numpy())[ltc:])
    ax[0].set_title('Inlet Pressure', fontdict={'fontsize': 24})
    ax[0].hlines(y = d2m(mPAP_sim), xmin = mpa['time'].to_numpy()[ltc], xmax = mpa['time'].to_numpy()[-1], linewidth=1, color='b')
    ax[1].plot(lpa['time'].to_numpy()[ltc:], d2m(lpa['pressure_out'].to_numpy())[ltc:])
    ax[1].set_title('LPA Outlet Pressure', fontdict={'fontsize': 24})
    ax[2].plot(rpa['time'].to_numpy()[ltc:], d2m(rpa['pressure_out'].to_numpy())[ltc:])
    ax[2].set_title('RPA Outlet Pressure', fontdict={'fontsize': 24})
    
    for i in range(3):
        ax[i].tick_params(axis="x", labelsize=16) 
        ax[i].tick_params(axis = 'y', labelsize=16)
        ax[i].set_xlabel('time (s)', fontdict={'fontsize': 20})
        ax[i].set_ylabel('pressure (mmHg)', fontdict={'fontsize': 20})
        
    ax[3].plot(mpa['time'].to_numpy()[ltc:], mpa['flow_in'].to_numpy()[ltc:])
    ax[3].set_title('Inlet Flow', fontdict={'fontsize': 24})
    ax[4].plot(lpa['time'].to_numpy()[ltc:], lpa['flow_out'].to_numpy()[ltc:])
    ax[4].set_title('LPA Outlet Flow', fontdict={'fontsize': 24})
    ax[5].plot(rpa['time'].to_numpy()[ltc:], rpa['flow_out'].to_numpy()[ltc:])
    ax[5].set_title('RPA Outlet Flow', fontdict={'fontsize': 24})
    
    for i in range(3, 6):
        ax[i].tick_params(axis="x", labelsize=16) 
        ax[i].tick_params(axis = 'y', labelsize=16)
        ax[i].set_xlabel('time (s)', fontdict={'fontsize': 20})
        ax[i].set_ylabel('flow (ml/s)', fontdict={'fontsize': 20})
    
    
    factor1 = 0
    if targets['maxPAP'][0] == targets['maxPAP'][1]:
        factor1 = 0.5
    ax[0].fill_between(x = mpa['time'].to_numpy()[ltc:], y1 = targets['maxPAP'][0]-factor1, y2 = targets['maxPAP'][1]+factor1, color = 'r', alpha = .3, label = f"Target Max PAP: ({targets['maxPAP'][0]}, {targets['maxPAP'][1]})")
    factor2 = 0
    if targets['minPAP'][0] == targets['minPAP'][1]:
        factor2 = 0.5
    ax[0].fill_between(x = mpa['time'].to_numpy()[ltc:], y1 = targets['minPAP'][0]-factor2, y2 = targets['minPAP'][1]+factor2, color = 'g', alpha = .3, label = f"Target Min PAP: ({targets['minPAP'][0]}, {targets['minPAP'][1]})")
    factor3 = 0
    if targets['mPAP'][0] == targets['mPAP'][1]:
        factor3 = 0.5
    ax[0].fill_between(x = mpa['time'].to_numpy()[ltc:], y1 = targets['mPAP'][0]-factor3, y2 = targets['mPAP'][1]+factor3, color = 'b', alpha = .3, label = f"Target Avg PAP: ({targets['mPAP'][0]}, {targets['mPAP'][1]})")
    
    ax[0].legend(fontsize = 16, loc = 'upper left', framealpha = .5)

    fig.savefig(str(save_dir / 'tuning_waveforms.png'))
    
##########
# Params #
##########
    
def load_tuning_params(TM: Manager):
    
    params = TuneParams()
    P = TM['tune_params']
    
    # loads params
    params.cap_wedge_pressure = m2d(P['PCWP'])
    params.minPAP_meas = [m2d(P["minPAP"][0]), m2d(P["minPAP"][1])]
    params.maxPAP_meas = [m2d(P["maxPAP"][0]), m2d(P["maxPAP"][1])]
    params.rpa_flow_split = P["rpa_split"]
    params.mPAP_meas = [m2d(P["mPAP"][0]), m2d(P["mPAP"][1])] 
    params.r_LPA = P["R_LPA"]
    params.r_RPA = P["R_RPA"]
    
    return params

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Tunes a model if necessary')
    
    parser.add_argument('-i', dest = 'config', help = 'Config.yaml file')
    parser.add_argument('--f', dest = 'force', action = 'store_true', default = False, help = 'Whether to restart tuning, even if tuning was already done once.')
    parser.add_argument('--s', dest = 'sensitivity_test', action = 'store_true', default = False, help = 'flag to run sensitivity tests or not')
    
    args = parser.parse_args()
    
    # manager
    TM = Manager(args.config)
    
    #! To be modified
    if args.sensitivity_test:
        raise NotImplementedError("Sensitivity Tests not implemented. Please remove the --s flag.")
    
    # check tune
    if not TM['options']['tune']:
        print('Tune option was set to false.')
        exit(1)
        
    # Create tuning directory
    tuning_dir = Path(TM['workspace']['lpn_dir']) / 'tuning'
    if args.force and tuning_dir.exists():
        shutil.rmtree(str(tuning_dir))

    try:
        tuning_dir.mkdir(exist_ok=False)
    except FileExistsError:
        print("Tuning dir already exists. Use --f flag to retune forcefully.")
        exit(1)
        
        
    print("Tuning Model " + TM['metadata']['model_name'] + "...")
    
    # load the main LPN
    main_lpn = LPN.from_file(str(TM['workspace']['lpn']))
    
    # add arguments to params
    params = load_tuning_params(TM)
    
    # set up tuning lpn
    tuning_lpn = construct_tuning_lpn(params, main_lpn)
    
    # run optimizer
    tune(TM, main_lpn, tuning_lpn, params, tuning_dir)
    
    # save a copy of the base lpn to main
    base_lpn_path = Path(TM['workspace']['root']) / 'base_lpn.in'
    TM.register('base_lpn', str(base_lpn_path), depth = ['workspace'])
    main_lpn.write_lpn_file(str(base_lpn_path))
