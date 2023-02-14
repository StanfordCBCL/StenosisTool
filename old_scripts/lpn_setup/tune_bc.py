# File: tune_bc.py
# File Created: Monday, 31st October 2022 8:46:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 1st February 2023 12:53:42 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Tunes Boundary Conditions for a 0D model using a simplified tuning model.
#! Sensitivity Tests to be implemented & clean up code



from sgt.core.solver import SolverResults, Solver0Dcpp
from sgt.utils.parser import ToolParser
from sgt.core.manager import TuningManager
from sgt.core.flow import Inflow
from sgt.core.bc import BoundaryConditions
from sgt.core.lpn import LPN
from sgt.utils.misc import m2d, d2m
from sgt.utils.io import write_json

from pathlib import Path
import numpy as np
import shutil
from copy import deepcopy
from scipy import optimize
from scipy.interpolate import interp1d
from tqdm import tqdm
from functools import partial
import matplotlib.pyplot as plt


##########
# Params #
##########


class TuneParams():
    
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
    
    ## vessels
    mpa_branch = main_lpn.get_mpa()
    
    ## mpa_junctions
    r = 0
    c = 0
    l = 0
    for vess in mpa_branch.vessel_info:
        r += vess['zero_d_element_values']['R_poiseuille']
        # compute sten coeff as R_expansion
        r += abs(main_lpn.inflow.mean_inflow) * vess['zero_d_element_values']['stenosis_coefficient']
        c += 1/vess['zero_d_element_values']['C']
        l += vess['zero_d_element_values']['L']
        
    c = 1/c

    mpa = [{
            "boundary_conditions": {
                "inlet": "INFLOW"
            },
            "vessel_id": 0,
            "vessel_length": 10.0,
            "vessel_name": "branch_mpa",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": r,
                "C": c,
                "L": l,
            }
        }
           ]
    
    
    
    lpa_res, rpa_res = compute_lpa_rpa_resistances(main_lpn)
    lpa =  [{
            "vessel_id": 1,
            "vessel_length": 10.0,
            "vessel_name": "branch_lpa_tree",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille":  lpa_res,
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
                "R_poiseuille": rpa_res,
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
    
    # check that the results are within reasonable bounds
    approximate_pdrop(main_lpn, None)

    
    return tuning_lpn

def approximate_pdrop(main_lpn: LPN, results: SolverResults = None):
    ''' approximates pressure drop accross LPN by computing total resistance in tree as one element'''
    tree = main_lpn.get_branch_tree()
    tree_val = compute_total_vals(tree, results)
    print('Approximate Pressure Drop from R values alone = ~', d2m(main_lpn.inflow.mean_inflow * tree_val[0]))

def compute_lpa_rpa_resistances(main_lpn: LPN, results: SolverResults = None):
    ''' computes LPA and RPA resistances from each branch, using resistances in series and parallel
    '''
    mpa = main_lpn.get_vessel_tree()
    while len(mpa.children) == 1:
        mpa = mpa.children[0]
        
    assert len(mpa.children) == 2, 'MPA should only branch into 2 sections, LPA and RPA'
        
    # determine which side child is.
    if mpa.children[0].side == 'lpa':
        lpa = mpa.children[0]
        rpa = mpa.children[1]
    else:
        lpa = mpa.children[1]
        rpa = mpa.children[0]
    
    # Compute RCL for each side.
    lpa_val = compute_total_vals(lpa, results)
    rpa_val = compute_total_vals(rpa, results)

    return lpa_val[0], rpa_val[0]

def compute_lpa_rpa_resistances2(main_lpn: LPN, results: SolverResults = None):
    ''' computes LPA and RPA resistances from each using mean pressure diff and flow.
    '''
    mpa = main_lpn.get_vessel_tree()
    while len(mpa.children) == 1:
        mpa = mpa.children[0]
        
    assert len(mpa.children) == 2, 'MPA should only branch into 2 sections, LPA and RPA'
        
    # determine which side child is.
    if mpa.children[0].side == 'lpa':
        lpa = mpa.children[0]
        rpa = mpa.children[1]
    elif mpa.children[0].side == 'rpa':
        lpa = mpa.children[1]
        rpa = mpa.children[0]

    # compute LPA approx
    mPAP_lpa = results.get_avg_val(lpa.vessel_info[0]['vessel_name'], 'pressure_in')
    Q_lpa = results.get_avg_val(lpa.vessel_info[0]['vessel_name'], 'flow_in')
    total_dP = []
    for node in main_lpn.tree_bfs_iterator(lpa):
        # boundary condition
        if "boundary_conditions" in node.vessel_info[-1]:
            tmp_mPAP_out = results.get_avg_val(node.vessel_info[-1]['vessel_name'], 'pressure_out')
            total_dP.append(mPAP_lpa - tmp_mPAP_out)
    print("LPA mean Inlet Pressure:", d2m(mPAP_lpa), "; LPA Average Pressure Drop:", d2m(np.array(total_dP).mean()))
    lpa_res = np.array(total_dP).mean() / Q_lpa
    
    # compute RPA approx
    mPAP_rpa = results.get_avg_val(rpa.vessel_info[0]['vessel_name'], 'pressure_in')
    Q_rpa = results.get_avg_val(rpa.vessel_info[0]['vessel_name'], 'flow_in')
    total_dP = []
    for node in main_lpn.tree_bfs_iterator(rpa):
        # boundary condition
        if "boundary_conditions" in node.vessel_info[-1]:
            tmp_mPAP_out = results.get_avg_val(node.vessel_info[-1]['vessel_name'], 'pressure_out')
            total_dP.append(mPAP_rpa - tmp_mPAP_out)
    print("RPA mean Inlet Pressure:", d2m(mPAP_rpa), "; RPA Average Pressure Drop:", d2m(np.array(total_dP).mean()))
    rpa_res = np.array(total_dP).mean() / Q_rpa
    
    return lpa_res, rpa_res


def compute_total_vals(cur_node: LPN.BranchNode, results: SolverResults = None):
    ''' computes total RCL coefficients of a tree recursively 
    Includes stenosis coefficients if results is not None
    '''
    if not cur_node.children:
        
        r = 0
        c = 0
        l = 0
        for vess in cur_node.vessel_info:
            
            # resistance factor
            r += vess['zero_d_element_values']['R_poiseuille']
            if results:
                # R expansion factor
                r += vess['zero_d_element_values']['stenosis_coefficient'] * abs(results.get_avg_val(vess['vessel_name'], 'flow_in'))
            
            # CL
            c += 1/vess['zero_d_element_values']['C']
            l += vess['zero_d_element_values']['L']
        
        c = 1/c
        return r, c, l
    
    r = 0
    c = 0
    l = 0
    # compute RCL in parallel for children.
    for child_node in cur_node.children:
        r_child, c_child, l_child = compute_total_vals(child_node, results)
        if r_child != 0:
            r += 1/r_child
        if c_child != 0:
            c += c_child
        if l_child != 0:
            l += 1/l_child
    r = 1/r
    l = 1/l
    c = 1/c
    
    for vess in cur_node.vessel_info:
            
        # resistance factor
        r += vess['zero_d_element_values']['R_poiseuille']
        if results:
            # R expansion factor
            r += vess['zero_d_element_values']['stenosis_coefficient'] * abs(results.get_avg_val(vess['vessel_name'], 'flow_in'))
        
        # CL
        c += 1/vess['zero_d_element_values']['C']
        l += vess['zero_d_element_values']['L']
    
    c = 1/c
    
    return r, c, l

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
    
    mPAP_loss, qRPA_loss, maxPAP_loss, minPAP_loss= loss_function(results = rez,
                         tune_params = params,
                         inflow = main_lpn.inflow,
                         )
    
    loss = mPAP_loss + qRPA_loss + maxPAP_loss + minPAP_loss
    
    
    params.iter += 1
    print(f"{params.iter:^15}|{mPAP_loss:^15.5f}|{maxPAP_loss:^15.5f}|{minPAP_loss:^15.5f}|{qRPA_loss:^15.5f}|{loss:^15.5f}")
    
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

def split_rcrs(TM: TuningManager, x, PCWP):
    ''' splits rcr driver '''
    capinfo = load_area_file(str(TM.capinfo))
    del capinfo[TM.inlet]
    
    rcrs = split_bc(capinfo, x)
    
    bcs = BoundaryConditions()
    
    for name, vals in rcrs.items():
        bcs.add_rcr(face_name=name, Rp = vals['Rp'], C = vals['C'], Rd = vals['Rd'], Pd = PCWP)
    return bcs


def add_rcrs(main_lpn: LPN, bcs: BoundaryConditions):
    ''' Adds appropriate rcrs to the main lpn
    '''
    
    bc_inv_map = bcs.get_bc_map()
    
    for bc in main_lpn.bc:
        if bc['bc_type'] == 'RCR':
            bc_values = bc_inv_map[main_lpn.bc_map[bc['bc_name']]]
            bc['bc_values']['Rp'] = bc_values['Rp']
            bc['bc_values']['Rd'] = bc_values['Rd']
            bc['bc_values']['C'] = bc_values['C']
            bc['bc_values']['Pd'] = bc_values['Pd']
    
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

#####################
# Update Tuning LPN #
#####################

def update_tuning_lpn(tuning_lpn: LPN, main_lpn: LPN, main_results: SolverResults):
    ''' updates the tuning param with new resistances include stenosis coefficients'''
    #approximate_pdrop(main_lpn, main_results)

    lpa_res, rpa_res = compute_lpa_rpa_resistances2(main_lpn, main_results)
    
    tuning_lpn._construct_vessel_map()
    lpa = tuning_lpn.get_vessel(1)
    lpa['zero_d_element_values']['R_poiseuille'] =  lpa_res

    rpa = tuning_lpn.get_vessel(4)
    rpa['zero_d_element_values']['R_poiseuille'] =  rpa_res
    
###############
# Tune Driver #
###############

def tune(TM: TuningManager, main_lpn: LPN, tuning_lpn: LPN, params: TuneParams, store_intermediate = False):
    ''' Tuning Steps
    '''
    # setup initial conditions
    x0 = get_initial_cond(params, main_lpn, tuning_lpn)
    modify_params(tuning_lpn, x0)
    tuning_lpn.write_lpn_file(TM.tuning_lpn)
    
    # split rcrs
    x_new = np.array([.2*x0[1],x0[0], .8*x0[1],.2*x0[3],x0[2], .8*x0[3] ])
    bcs = split_rcrs(TM, x_new, params.cap_wedge_pressure)
    add_rcrs(main_lpn, bcs)

    # run a simulation on the main lpn
    main_solver = Solver0Dcpp(main_lpn, use_steady=True, last_cycle_only=True, debug = True)
    main_results = main_solver.run_sim()
        
    # update Tuning LPN accordingly
    update_tuning_lpn(tuning_lpn, main_lpn, main_results)
    
    
    # bounds
    bounds = optimize.Bounds([0,0,0,0], [ np.inf, np.inf, np.inf, np.inf], keep_feasible=True)
    
    # run it max epoch # of times or till stable.
    diff_tol = .01 # 1 %
    max_epochs = 10
    epoch = 1
    diff = 1
    prev = None
    cur = None
    while epoch < max_epochs:
        # run optimizer
        params.iter = 0 
        print(f"{'Iteration':^15}|{'mPAP Loss':^15}|{'maxPAP Loss':^15}|{'minPAP Loss':^15}|{'qRPA Loss':^15}|{'Total Loss':^15}")
        print("-" * 15 * 7 + "-" * 6)
        results = optimize.minimize(fun = opt_function,
                                    x0 = x0,
                                    args = (main_lpn,
                                            tuning_lpn, params),
                                    method='Nelder-Mead',
                                    bounds=bounds,
                                    options = {'disp': True,
                                               'maxiter': 200})
            
        # set prev x0 as start point
        x0 = results.x
        # split using 1:4 ratio
        x_new = np.array([.2*x0[1],x0[0], .8*x0[1],.2*x0[3],x0[2], .8*x0[3] ])
        
        
        # track difference.
        prev = cur
        cur = x0
        print("Curr", cur)
        if prev is not None and cur is not None:
            print("Prev:", prev)
            diff = np.nan_to_num(abs(prev - cur) / prev, nan =0 ).sum()
            print("Difference from previous Epoch: ", diff)
            if diff <= diff_tol:
                break
            

        # store intermediate x values if requested
        if store_intermediate:
            epoch_dir = TM.intermediate_dir / f'epoch{epoch}'
            epoch_dir.mkdir(exist_ok = True)
            # save x
            np.save(str(epoch_dir / 'results.npy'), results.x, allow_pickle=True)
            # save tuning file
            modify_params(tuning_lpn, results.x)
            tuning_lpn.write_lpn_file(epoch_dir / 'tuner.in')
            # validate results
            validate_results(params,
                             main_lpn,
                             tuning_lpn,
                             x_new,
                             epoch_dir,
                             TM.tune_params
                             )
                    
        
        # split rcrs
        bcs = split_rcrs(TM, x_new, params.cap_wedge_pressure)
        add_rcrs(main_lpn, bcs)
    
        # run a simulation on the main lpn
        main_solver = Solver0Dcpp(main_lpn, use_steady=True, last_cycle_only=True, debug = True)
        main_results = main_solver.run_sim()
            
        # update Tuning LPN accordingly
        update_tuning_lpn(tuning_lpn, main_lpn, main_results)
        epoch += 1
   
        
    
    # store final
    np.save(str(TM.tuning_dir / 'results.npy'), results.x, allow_pickle=True)
    # save tuning file
    modify_params(tuning_lpn, results.x)
    tuning_lpn.write_lpn_file(TM.tuning_lpn)
    
    x0 = results.x
    x_new = np.array([.2*x0[1],x0[0], .8*x0[1],.2*x0[3],x0[2], .8*x0[3] ])
    # split rcrs
    bc = split_rcrs(TM, x_new, params.cap_wedge_pressure)
    add_rcrs(main_lpn, bc)
    
    # run a simulation on the main lpn
    main_solver = Solver0Dcpp(main_lpn, use_steady=True, last_cycle_only=False, debug = True)
    main_results = main_solver.run_sim()
    main_results.validate_results(main_lpn, outfile=str(TM.tuning_dir / "full_inlet_waveform.png"), targets = TM.tune_params)
    
    # validate results
    validate_results(params,
                        main_lpn,
                        tuning_lpn,
                        x_new,
                        TM.tuning_dir,
                        TM.tune_params)
    
    print("Optimization Complete.")
    
    return bc
    
    
    

def get_initial_cond(params: TuneParams, main_lpn: LPN, tuning_lpn: LPN):
    ''' initial cond with numpy array in form of [Rp_LPA, C_LPA, Rd_LPA, Rp_RPA, C_RPA, Rd_RPA]'''
    # determine PVR
    PVR = ((params.mPAP_meas[0] + params.mPAP_meas[1])/2 - params.cap_wedge_pressure)/ main_lpn.inflow.mean_inflow
    
    # determine MPA capacitance using capacitance in series
    mpa_c = tuning_lpn.vessel[0]['zero_d_element_values']['C']
    
    x0 = np.array([
                    mpa_c,
                    2 * PVR,
                    mpa_c,
                    2 * PVR]
                    )
    return x0

def modify_params(lpn: LPN, x ):
    ''' modifies the optimizer variables '''
    vess = lpn.vessel

    # LPA
    vess[2]['zero_d_element_values']['C'] = x[0]
    vess[3]['zero_d_element_values']['R_poiseuille'] = x[1]
    
    # RPA
    vess[5]['zero_d_element_values']['C'] = x[2]
    vess[6]['zero_d_element_values']['R_poiseuille'] = x[3]

def convert_to_dict(opt_results: optimize.OptimizeResult):
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

####################
# Sensitivity Test #
####################

# def sens_test(x, tune_params : TuneParams, inflow: Inflow0D, solver: Solver0D, sensitivity_dir):
#     mapping = {'Rp_LPA': 0, 'C_LPA': 1, 'Rd_LPA': 2, 'Rp_RPA':3, 'C_RPA':4, 'Rd_RPA':5}
    
#     for var_name, index in mapping.items():
#         print(f'\tRunning {var_name} sensitivity test...', end = '\t', flush = True)
        
#         fig, ax = plt.subplots(2, 2, figsize = (40, 30))
#         ax = ax.flatten()
#         ax[0].set_xlabel(f'pct change {var_name}', fontdict={'fontsize': 20})
#         ax[1].set_xlabel(f'pct change {var_name}', fontdict={'fontsize': 20})
#         ax[0].set_ylabel(f'mPAP', fontdict={'fontsize': 20})
#         ax[1].set_ylabel(f'Q_RPA_avg', fontdict={'fontsize': 20})
#         ax[0].set_title('mPAP change', fontdict={'fontsize': 24})
#         ax[1].set_title('Q_RPA change', fontdict={'fontsize': 24})
#         ax[2].set_xlabel(f'pct_change {var_name}', fontdict={'fontsize': 20})
#         ax[2].set_ylabel(f'MSE', fontdict={'fontsize': 20})
#         ax[2].set_title('MSE', fontdict={'fontsize': 24})
#         ax[3].set_xlabel(f'pct_change {var_name}', fontdict = {'fontsize':20})
#         ax[3].set_ylabel(f'maxPAP', fontdict={'fontsize': 20})
#         ax[3].set_title('maxPAP', fontdict={'fontsize': 24})
        
#         for i in range(3):
#             ax[i].tick_params(axis="x", labelsize=16) 
#             ax[i].tick_params(axis = 'y', labelsize=16)
        
#         mod = np.linspace(.8, 1.2, 40)
#         mpap = []
#         q_rpa = []
#         mses = []
#         maxpap = []
#         for pct in mod:
#             x_test = np.copy(x)
#             x_test[index] *= pct
            
#             rez = run_sim_wrapper(x_test, solver)
            
#             _, _, mse, _ , (mPAP_sim, Q_RPA_sim, maxPAP_sim, _, _,_) = loss_function(rez, tune_params, inflow)
#             mpap.append(mPAP_sim)
#             q_rpa.append(Q_RPA_sim)
#             mses.append(mse)
#             maxpap.append(maxPAP_sim)
            

        
#         ax[0].plot(mod - 1, mpap)
#         ax[1].plot(mod - 1, q_rpa)
#         ax[2].plot(mod - 1, mses)
#         ax[3].plot(mod - 1, maxpap)
        
#         fig.savefig(os.path.join(sensitivity_dir, f'{var_name}_change.png'))
#         print('Done')
        
###############
# Main Driver #
###############

def main(TM: TuningManager, args):
    
    print("Tuning Model " + TM.model_name + "...")
    
    # load the main LPN
    main_lpn = LPN.from_file(str(TM.lpn))
    
    #? Implement Config here if necessary
    tuning_params = TuneParams()
    #add_config(tuning_params, args)
    add_config(tuning_params, TM)
    
    # set up tuning lpn
    tuning_lpn = construct_tuning_lpn(tuning_params, main_lpn)
    
    # run optimizer
    bc = tune(TM, main_lpn, tuning_lpn, tuning_params, args.store_intermediate)
    # write the BC
    bc.write_rcrt_file(str(TM.rcrt.parent))
    # write lpn dir
    main_lpn.write_lpn_file(TM.lpn)
    
    # set tune to false
    TM.config_add(['options', 'tune'], False)
    TM.write_config()
    
# def add_config(params: TuneParams, args):
    
#     if args.PCWP:
#         params.cap_wedge_pressure = m2d(args.PCWP)
#     if args.minPAP:
#         params.minPAP_meas = [m2d(args.minPAP), m2d(args.minPAP)]
#     if args.maxPAP:
#         params.maxPAP_meas = [m2d(args.maxPAP), m2d(args.maxPAP)]
#     if args.rpa_split:
#         params.rpa_flow_split = args.rpa_split
#     if args.mPAP:
#         params.mPAP_meas = [m2d(args.mPAP), m2d(args.mPAP)]
def add_config(params: TuneParams, TM: TuningManager):
    
    params.cap_wedge_pressure = m2d(TM.tune_params['PCWP'])
    params.minPAP_meas = [m2d(TM.tune_params["minPAP"][0]), m2d(TM.tune_params["minPAP"][1])]
    params.maxPAP_meas = [m2d(TM.tune_params["maxPAP"][0]), m2d(TM.tune_params["maxPAP"][1])]
    params.rpa_flow_split = TM.tune_params["rpa_split"]
    params.mPAP_meas = [m2d(TM.tune_params["mPAP"][0]), m2d(TM.tune_params["mPAP"][1])] 




if __name__ == '__main__':

    tool = ToolParser(desc = 'Tunes a model if necessary')
    
    # dev params
    tool.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Whether to run an optimization even if tuning is set to false.')
    tool.parser.add_argument('-i', dest = 'store_intermediate', action = 'store_true', default = False, help = 'save intermediate results')
    tool.parser.add_argument('-s', dest = 'sensitivity_test', action = 'store_true', default = False, help = 'flag to run sensitivity tests or not')
    
    tool.parser.add_argument('-PCWP', type = float, default = None, help = 'Capillary Wedge pressure (mmHg)')
    tool.parser.add_argument('-mPAP', type = float, default = None, help = 'Mean Pulmonary Arterial pressure (mmHg)')
    tool.parser.add_argument('-maxPAP', type = float, default = None, help = 'Systolic Pulmonary Arterial Pressure (mmHg)')
    tool.parser.add_argument('-minPAP', type = float, default = None, help = 'Diastolic Pulmonary Arterial Pressure (mmHg)')
    tool.parser.add_argument('-rpa_split', type = float, default = None, help = 'fraction of flow going to Right Pulmonary Artery (0-1)')
    
    args = tool.parse_args()
    
    
    TM = TuningManager(args.config)
    
    if args.sensitivity_test:
        raise NotImplementedError("Sensitivity Tests not implemented. Please remove the -s flag.")
    
    if not TM.tune:
        if args.force:
            # clear previous TM values
            print("Removing previous existing tuning results...", end = '\t', flush = True)
            shutil.rmtree(str(TM.tuning_dir))
            print("Done")
            # recreate dirs
            TM = TuningManager(args.config)
        else:
            raise ValueError('Tune was set to false. Use -f if you wish to tune anyways.')
    
    main(TM, args)
        
        
    
    