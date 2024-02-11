
import json
from collections import defaultdict, deque, OrderedDict
import numpy as np
from typing import Generator, Union
from pathlib import Path
from copy import deepcopy
from scipy.interpolate import interp1d

from svinterface.utils.io import write_json
from svinterface.core.bc import Inflow, RCR
from svinterface.core.polydata import Centerlines
from abc import ABC, abstractclassmethod


class FastLPN():
    """A fast handler for LPN with minimal vessel and junction modification
    """
    # config file params
    BC = 'boundary_conditions'
    VESS = 'vessels'
    SIM = 'simulation_parameters'
    JUNC = 'junctions'
    DESC = 'description'
    BC_MAP = 'bc_map'
    
    
    def __init__(self, lpn_data):
        self.lpn_data = lpn_data
        
    @classmethod
    def from_file(cls, lpn_file):
        '''  loads LPN from a file
        '''
        with open(lpn_file, "r") as lfile:
            return cls(lpn_data = json.load(lfile))
        
    def copy(self):
        """Deep copy of itself

        Returns:
            FastLPN: a copy of this LPN
        """
        return FastLPN(deepcopy(self.lpn_data))

    def get_junction(self, id):
        """Gets junction according to ordering

        Args:
            id (int): junction id
        """
        return self.lpn_data[self.JUNC][id]
        
    def get_vessel(self, id):
        """Retrieves vessel

        Args:
            id (int): vessel_id

        Returns:
            dict: vessel dictionary values to be modified.
        """
        return self.lpn_data[self.VESS][id]
    
    def change_vessel(self, vessel_id: int, R: float = None, C: float = None, L: float = None, S: float = None, mode: str = 'replace'):
        """Changing Vessel

        Args:
            vessel_id (int): junction id
            R (float, optional): R poiseuille resistance. Defaults to None.
            C (float, optional): capacitance. Defaults to None.
            L (float, optional): inductance. Defaults to None.
            S (float, optional): stenosis coefficient. Defaults to None.
            mode (str, optional): replace or add. Defaults to replace
        """
        vess = self.get_vessel(vessel_id)
        vals = vess['zero_d_element_values']
        if R is not None:
            vals['R_poiseuille'] = R if mode == 'replace' else R + vals['R_poiseuille']
        if C is not None:
            vals['C'] = C if mode == 'replace' else C + vals['C']
        if L is not None:
            vals['L'] = L if mode == 'replace' else L + vals['L']
        if S is not None:
            vals['stenosis_coefficient'] = S if mode == 'replace' else S + vals['stenosis_coefficient']
    
    def change_junction_outlet(self, junction_id_or_name: int, which: int, R: float = None, C: float = None, L: float = None, S: float = None, mode: str = 'replace'):
        """Changing Junction Outlets

        Args:
            junction_id (int): junction id
            which (int): which outlet of the junction to change
            R (float, optional): R poiseuille resistance. Defaults to None.
            C (float, optional): capacitance. Defaults to None.
            L (float, optional): inductance. Defaults to None.
            S (float, optional): stenosis coefficient. Defaults to None.
            mode (str, optional): replace or add. Defaults to replace
        """
        junc = self.get_junction(junction_id_or_name)
        assert which < len(junc['outlet_vessels']), f"Selected outlet {which} is out of bounds"
        vals = junc['junction_values']
        if R is not None:
            vals['R_poiseuille'][which] = R if mode == 'replace' else R + vals['R_poiseuille'][which]
        if C is not None:
            vals['C'][which] = C if mode == 'replace' else C + vals['C'][which]
        if L is not None:
            vals['L'][which] = L if mode == 'replace' else L + vals['L'][which]
        if S is not None:
            vals['stenosis_coefficient'][which] = S if mode == 'replace' else S + vals['stenosis_coefficient'][which]
    
    def occlude_vessel(self, vessel_id: int, occlusion: float):
        ''' creates an occlusion in a vessel
        '''
        vess = self.get_vessel(vessel_id)
        remaining_area_pct = (1 - occlusion)
        # modify elements
        vess['zero_d_element_values']['R_poiseuille'] /= (remaining_area_pct**2)
        vess['zero_d_element_values']['C'] *= remaining_area_pct
        vess['zero_d_element_values']['L'] /= remaining_area_pct
        vess['zero_d_element_values']['stenosis_coefficient'] /= (remaining_area_pct**2)
        return vess
        
    
    def repair_vessel(self, vessel_id: int, area_increase: float):
        ''' repairs a vessel by an area_increase amount.
        '''
        vess = self.get_vessel(vessel_id)
        area_increase += 1
        # modify elements
        vess['zero_d_element_values']['R_poiseuille'] /= (area_increase**2)
        vess['zero_d_element_values']['C'] *= area_increase
        vess['zero_d_element_values']['L'] /= area_increase
        vess['zero_d_element_values']['stenosis_coefficient'] /= (area_increase**2)
        return vess

class OriginalLPN():
    """LPN for the original model creation.
    """
    BC = 'boundary_conditions'
    VESS = 'vessels'
    SIM = 'simulation_parameters'
    JUNC = 'junctions'
    DESC = 'description'
    
    def __init__(self):
        # main data
        self.lpn_file = None
        self.lpn_data = None
        
        # abstracted data
        self.inflow = None
        self.bc_data = None
        self.vessel_map = None
        self.junction_map = None
        
        
    ############
    # Init Ops #
    ############

    @classmethod
    def from_file(cls, lpn_file):
        '''  loads LPN from a file
        '''
        lpn = cls()
        lpn.read_lpn_file(lpn_file)
        return lpn
        
    
    @classmethod
    def from_dict(cls, solver_dict):
        ''' loads LPN from a dict
        '''
        lpn = cls()
        lpn.lpn_data = solver_dict
        lpn._update_lpn_data()
        return lpn
    
    def setup_empty_lpn(self):
        ''' sets up a new lpn
        '''
        self.lpn_data = {self.BC: [], 
                        self.VESS: [],
                        self.SIM: {},
                        self.JUNC: [],
                        self.DESC: {}}
        self._update_lpn_data()
    
    ########################
    # LPN Creation Methods #
    ########################
    
    def get_original_lpn(self):
        ''' returns a deep copy of the originalLPN'''
        return OriginalLPN.from_dict(deepcopy(self.lpn_data))
    
    def get_fast_lpn(self):
        """Retrieves a FastLPN version of the lpn.
        """
        return FastLPN(deepcopy(self.lpn_data))
    
    def get_full_lpn(self):
        """Initializes a full LPN
        """
        return LPN.from_dict(deepcopy(self.lpn_data))
    
    def to_cpp(self):
        """Converts to cpp by converting to BloodVesselJunctions
        """
        self.normal_to_bvj()
        
    def to_python(self):
        ''' converts BVJ to Normal Junctions'''
        for junc in self.junctions:
            if junc['junction_type'] == 'BloodVesselJunction':
                junc['junction_type'] = 'NORMAL_JUNCTION'
        
    def normal_to_bvj(self):
        """ Converts Normal junctions to BVJ junctions"""
        for junc in self.junctions:
            if junc['junction_type'] == 'NORMAL_JUNCTION':
                junc['junction_type'] = 'BloodVesselJunction'
                junc['junction_values'] = {"R_poiseuille": [0 for i in range(len(junc['outlet_vessels']))],
                                            "C": [0 for i in range(len(junc['outlet_vessels']))],
                                            "L": [0 for i in range(len(junc['outlet_vessels']))],
                                            "stenosis_coefficient": [0 for i in range(len(junc['outlet_vessels']))]}
    
    ##########
    # IO Ops #
    ##########
        
    def write_lpn_file(self, lpn_file: Path):
        ''' writes a dict into the lpn file
        '''
        write_json(lpn_file, self.lpn_data)
    
    def read_lpn_file(self, lpn_file):
        ''' reads the solver file into a dict
        '''
        self.lpn_file = lpn_file
        with open(lpn_file, 'r') as sfile:
            self.lpn_data = json.load(sfile)
        self._update_lpn_data()
    
    ###########################
    # Private Background Meth #
    ###########################
    
    def _update_lpn_data(self):
        ''' update attributes if solver data changes
        '''
        # compute inflow if it exists. Otherwise use a dummy.
        self.bc_data = OrderedDict()
        for bc in self.bc:
            if bc['bc_type'] == 'FLOW':
                self.inflow = Inflow(inflow_arr = np.array(list(zip(bc['bc_values']['t'], bc['bc_values']['Q']))))
            else:
                self.bc_data[bc['bc_name']] = bc
                
        if self.inflow is None:
            self.inflow = Inflow(inflow_arr = np.array([(0,0),(1,0)]), smooth = False)
            
        # construct vessel maps and junction matrix
        self._construct_vessel_maps()

    def _construct_vessel_maps(self):
        ''' Vessel & Junction Map to avoid Linear searches
        '''
        self.vessel_map = {}
        self.junction_map = {}
        for vess in self.vessel:
            self.vessel_map[vess['vessel_name']] = vess
        for junction in self.junctions:    
            self.junction_map[junction['junction_name']] = junction
    
    #################
    # Basic Methods #
    #################
    
    def get_junction(self, name_or_id):
        """Gets junction by name or id
        """
        if type(name_or_id) == str:
            return self.junction_map[name_or_id]
        elif type(name_or_id) == int:
            return self.junctions[name_or_id]
        
    def get_vessel(self, name_or_id):
        """Retrieves vessel by name or id
        """
        if type(name_or_id) == str:
            return self.vessel_map[name_or_id]
        elif type(name_or_id) == int:
            return self.vessel[name_or_id]
        
    def get_outlet_bc(self, bc_name):
        """Retrieves outlet BC by name
        """
        return self.bc_data[bc_name]
    
    def get_inlet_vessel(self):
        ''' identify the inlet vessel
        '''
        for vess in self.vessel:
            if 'boundary_conditions' in vess:
                if 'inlet' in vess['boundary_conditions']:
                    return vess['vessel_id'], vess['vessel_name']
        print('Error: an inlet bc does not exist.')
        return None
    
    def change_vessel(self, vessel_id_or_name: int, R: float = None, C: float = None, L: float = None, S: float = None, mode: str = 'replace'):
        """Changing Vessel

        Args:
            vessel_id (int): junction id
            R (float, optional): R poiseuille resistance. Defaults to None.
            C (float, optional): capacitance. Defaults to None.
            L (float, optional): inductance. Defaults to None.
            S (float, optional): stenosis coefficient. Defaults to None.
        """
        vess = self.get_vessel(vessel_id_or_name)
        vals = vess['zero_d_element_values']
        if R is not None:
            vals['R_poiseuille'] = R if mode == 'replace' else R + vals['R_poiseuille']
        if C is not None:
            vals['C'] = C if mode == 'replace' else C + vals['C']
        if L is not None:
            vals['L'] = L if mode == 'replace' else L + vals['L']
        if S is not None:
            vals['stenosis_coefficient'] = S if mode == 'replace' else S + vals['stenosis_coefficient']
    
    def change_junction_outlet(self, junction_id_or_name: int, which: int, R: float = None, C: float = None, L: float = None, S: float = None, mode: str = 'replace'):
        """Changing Junction Outlets

        Args:
            junction_id (int): junction id
            which (int): which outlet of the junction to change
            R (float, optional): R poiseuille resistance. Defaults to None.
            C (float, optional): capacitance. Defaults to None.
            L (float, optional): inductance. Defaults to None.
            S (float, optional): stenosis coefficient. Defaults to None.
        """
        junc = self.get_junction(junction_id_or_name)
        assert which < len(junc['outlet_vessels']), f"Selected outlet {which} is out of bounds"
        assert junc['junction_type'] == "BloodVesselJunction", f"Selected Junction {junction_id_or_name} is not a BloodVesselJunction."
        vals = junc['junction_values']
        if R is not None:
            vals['R_poiseuille'][which] = R if mode == 'replace' else R + vals['R_poiseuille'][which]
        if C is not None:
            vals['C'][which] = C if mode == 'replace' else C + vals['C'][which]
        if L is not None:
            vals['L'][which] = L if mode == 'replace' else L + vals['L'][which]
        if S is not None:
            vals['stenosis_coefficient'][which] = S if mode == 'replace' else S + vals['stenosis_coefficient'][which]
    
    def occlude_vessel(self, vessel_id_or_name: int, occlusion: float):
        ''' creates an occlusion in a vessel
        '''
        vess = self.get_vessel(vessel_id_or_name)
        remaining_area_pct = (1 - occlusion)
        # modify elements
        vess['zero_d_element_values']['R_poiseuille'] /= (remaining_area_pct**2)
        vess['zero_d_element_values']['C'] *= remaining_area_pct
        vess['zero_d_element_values']['L'] /= remaining_area_pct
        vess['zero_d_element_values']['stenosis_coefficient'] /= (remaining_area_pct**2)
        return vess
        
    
    def repair_vessel(self, vessel_id_or_name: int, area_increase: float):
        ''' repairs a vessel by an area_increase amount.
        '''
        vess = self.get_vessel(vessel_id_or_name)
        area_increase += 1
        # modify elements
        vess['zero_d_element_values']['R_poiseuille'] /= (area_increase**2)
        vess['zero_d_element_values']['C'] *= area_increase
        vess['zero_d_element_values']['L'] /= area_increase
        vess['zero_d_element_values']['stenosis_coefficient'] /= (area_increase**2)
        return vess
    
    def get_vessel_radius(self, vessel_id_or_name: int):
        ''' reverse engineer radius from resistance
        '''
        vess = self.get_vessel(vessel_id_or_name)
        resistance = vess['zero_d_element_values']['R_poiseuille'] 
        radius = (8 * self.simulation_params['viscosity'] * vess['vessel_length'] / (np.pi + resistance)) ** (1/4)
        return radius
        
    def num_vessels(self):
        return len(self.vessel)
    
    def num_junctions(self):
        return len(self.junctions)
    
    ##############
    # Properties #
    ##############
    
    @property
    def vessel(self):
        return self.lpn_data[self.VESS]

    @vessel.setter
    def vessel(self, val):
        self.lpn_data[self.VESS] = val
        self._update_lpn_data()
        
    @property
    def bc(self):
        return self.lpn_data[self.BC]

    @bc.setter
    def bc(self, val):
        self.lpn_data[self.BC] = val
        self._update_lpn_data()
    
    @property
    def simulation_params(self):
        return self.lpn_data[self.SIM]

    @simulation_params.setter
    def simulation_params(self, val):
        self.lpn_data[self.SIM] = val
        self._update_lpn_data()

    @property
    def junctions(self):
        return self.lpn_data[self.JUNC]

    @junctions.setter
    def junctions(self, val):
        self.lpn_data[self.JUNC] = val
        self._update_lpn_data()
        
    @property
    def description(self):
        return self.lpn_data[self.DESC]

    @description.setter
    def description(self, val):
        self.lpn_data[self.DESC] = val
        self._update_lpn_data()
    

class LPN(OriginalLPN):
    ''' Full LPN File handler with versatile option to get data. Slowest & not meant for speed computations '''

    FLAGS = 'flags'
    
    class Node(ABC):
        ''' base node for trees
        '''
        def __init__(self, ids, vessel_info):
            self.ids = ids
            self.vessel_info = vessel_info
            self.parent = None
            self.children = []
            
        def last_vessel(self):
            ''' get the last vessel
            '''
            return self.ids[-1], self.vessel_info[-1]
        
        def set_metadata(self, key, value):
            ''' sets metadata for all vessels at once'''
            for vess in self.vessel_info:
                vess[key] = value
        
        def __len__(self):
            return len(self.ids)
        
        @abstractclassmethod
        def type(self):
            pass
        
        @abstractclassmethod
        def id(self):
            pass
        
        @abstractclassmethod
        def length():
            pass

    class BranchNode(Node):
        ''' node for branch tree
        '''
        def __init__(self, vess_id: list, vessel_info: list):
            super().__init__(ids=[vess_id], vessel_info=[vessel_info])
        
        @property
        def type(self):
            return "branch"
        
        @property
        def id(self):
            """Branch ID"""
            return int(self.vessel_info[-1]['vessel_name'].split('_')[0][6:])
        
        def length(self):
            ''' compute total length of the branch
            '''
            length = 0
            for vessel in self.vessel_info:
                length += vessel['vessel_length']
            return length
    
    class JunctionNode(Node):
        
        def __init__(self, junction_name, vessel_info):
            super().__init__(ids = [junction_name], vessel_info=[vessel_info])

        @property
        def type(self):
            return "junction"
        
        @property
        def id(self):
            """Junction ID"""
            return int(self.vessel_info[-1]['junction_name'][1:])
    
        def length(self):
            ''' all lengths of the junction
            '''
            return self.vessel_info[-1]['lengths']
        
    # flags we use for additional data
    FLAGS_PRESET = {"rcrt_map": False,
                    "sides": False,
                    "gid": False}
    
    
    def __init__(self):
        super().__init__()
    
    def __getattr__(self, item):
        return None
    
    def _update_lpn_data(self):
        if self.FLAGS not in self.lpn_data:
            self.lpn_data[self.FLAGS]=self.FLAGS_PRESET
        super()._update_lpn_data()
        
    @property
    def flags(self):
        return self.lpn_data[self.FLAGS]
    
    def update(self):
        self.write_lpn_file(self.lpn_file)
    
    ###################
    # Additional Data #
    ###################
    
    def add_rcrt_map(self, outlet_faces: list, overwrite = True):
        ''' must be in same order as bc that was written to rcrt.dat file (same order as outlet_face_names.dat). Automatically flushed.
        '''
        
        if self.flags['rcrt_map'] and not overwrite:
            raise ValueError("RCRTs were already mapped once.")
        
        # for each bc right now,
        for idx, face in enumerate(self.bc_data):
            self.bc_data[face]['face_name'] = outlet_faces[idx]
        self.flags['rcrt_map'] = True
        
    def update_rcrs(self, rcrs: RCR):
        """ Updates RCRs from an rcr class"""
        if not self.flags['rcrt_map']:
            raise ValueError("RCR_map must be set before rcrs can be used")
        
        for bc in self.bc:
            if 'face_name' in bc:
                updated_bc = rcrs.bc_list[bc['face_name']]
                bc['bc_values']['Rp'] = updated_bc['Rp']
                bc['bc_values']['C'] = updated_bc['C']
                bc['bc_values']['Rd'] = updated_bc['Rd']
                bc['bc_values']['Pd'] = updated_bc['Pd']
                
        
            
    def det_lpa_rpa(self, head_node: Node, overwrite = True):
        ''' Determine whether a node is an LPA or an RPA
        '''
        # check if run already
        if self.flags['sides'] and not overwrite:
            print("Sides has already been run. To overwrite, use overwrite = True")
            return
        
        # needs to have rcrt_map
        if not self.flags['rcrt_map']:
            print("rcrt_map flag must be set to true: try running add_rcrt_map first.")
            return 
        
        # sides of PA
        head_node.set_metadata('side', 'MPA')
        
        junc1 = head_node.children[0]
        junc1.set_metadata('side', 'MPA')
        if len(junc1.children) != 2:
            print("MPA can only have 2 children, LPA and RPA. Please refine your mesh/centerlines.")
            return
        
        side1 = junc1.children[0]
        side2 = junc1.children[1]
        
        # iterate down until first child is reached
        cur_branch = side1
        while cur_branch.children:
            cur_branch = cur_branch.children[0]

        # Check if rcr contains 'lpa' or 'rpa' & fill appropriately
        rcr_name = self.bc_data[cur_branch.vessel_info[-1]['boundary_conditions']['outlet']]['face_name']
        if 'lpa' in rcr_name.lower():
            for node in self.tree_bfs_iterator(side1):
                node.set_metadata('side', 'LPA')
            for node in self.tree_bfs_iterator(side2):
                node.set_metadata('side', 'RPA')
        elif 'rpa' in rcr_name.lower():
            for node in self.tree_bfs_iterator(side1):
                node.set_metadata('side', 'RPA')
            for node in self.tree_bfs_iterator(side2):
                node.set_metadata('side', 'LPA')
        else:
            print('Outlet names must all contain "lpa" and "rpa", so a side could not be determined.')
            return 
        
        self.flags['sides'] = True
        

    def find_gids(self, cent: Centerlines, overwrite= True):
        """ Finds the global node ID's of where the junctions are on centerlines.
        #! Centerlines MUST be the same as the ones used to create the LPN, otherwise the ids will be off
        """
          # check if run already
        if self.flags['gid'] and not overwrite:
            print("Gid has already been run. To overwrite, use overwrite = True")
            return
        
        tree = self.get_tree()
        # get array data.
        br_id = cent.get_pointdata_array(cent.PointDataFields.BRANCHID)
        g_id = cent.get_pointdata_array(cent.PointDataFields.NODEID)
        paths = cent.get_pointdata_array(cent.PointDataFields.PATH)
        
        
        valid_junctions = np.zeros(g_id.shape[0], dtype=int) - 1
        valid_vessels = np.zeros(g_id.shape[0], dtype=int) - 1
        valid_caps = np.zeros(g_id.shape[0],dtype=int) - 1
        
        for node in self.tree_bfs_iterator(tree, allow="branch"):
            
            # branch node
            br = node.id
            br_gid = g_id[br_id == br]
            br_paths = paths[br_id == br]
        
            
            # initial vessel
            vess_in = 0
            vess_out = node.vessel_info[0]['vessel_length']
            gid_in = br_gid[np.argmin(abs(br_paths - vess_in))].item()
            gid_out = br_gid[np.argmin(abs(br_paths - vess_out))].item()
            # append gid in and out
            node.vessel_info[0]['gid'] = [gid_in, gid_out]
            valid_vessels[gid_in] = node.vessel_info[0]['vessel_id']
            valid_vessels[gid_out] = node.vessel_info[0]['vessel_id']
            
            
            
            for vess in node.vessel_info[1:]:
                # following vessels
                vess_in = vess_out
                vess_out += vess['vessel_length']
                # allow for slight inaccuracies
                gid_in = br_gid[np.argmin(abs(br_paths - vess_in))].item()
                gid_out = br_gid[np.argmin(abs(br_paths - vess_out))].item()
                # append gid in and out
                vess['gid'] = [gid_in, gid_out]
                valid_vessels[gid_in] = vess['vessel_id']
                valid_vessels[gid_out] = vess['vessel_id']

            # add corresponding vessel the cap belongs to.
            last_vess = node.vessel_info[-1]
            if 'boundary_conditions' in last_vess:
                if 'outlet' in last_vess['boundary_conditions']:
                    valid_caps[last_vess['gid'][-1]] = last_vess['vessel_id']
                
        # once all branches have been resolved
        for node in self.tree_bfs_iterator(tree, allow = "junction"):
            jc = node.id
            
            # upstream vessel's out gid
            gid_in = node.parent.vessel_info[-1]['gid'][-1]
            # downstream vessel's in gid
            gid_out = [c.vessel_info[0]['gid'][0] for c in node.children]
            
            # downstream vessels
            node.vessel_info[0]['gid'] = [gid_in, gid_out]
            
            # update valid junctions
            valid_junctions[gid_in] = node.id
            valid_junctions[gid_out] = node.id
        
        cent.add_pointdata(valid_junctions, "Junctions_0D")
        cent.add_pointdata(valid_caps, "Caps_0D")
        cent.add_pointdata(valid_vessels, "Vessels_0D")
        
        self.flags['gid'] = True
        return cent
    
    #####################
    # Tree Construction #
    #####################
    
    def _junction_matrix(self):
        ''' Junction Maxtrix to map inlet and outlet of junctions
        '''
        junc_mat = defaultdict(lambda: None)
        for junc in self.junctions:
            junc_mat[junc['inlet_vessels'][0]] = junc
        return junc_mat
    
    def get_tree(self) -> Node:
        ''' get the tree at that snapshot in the class
        '''
        
        # construct junction matrix
        junc_mat = self._junction_matrix()
        
        head_node = self.BranchNode(vess_id=0, 
                                    vessel_info=self.get_vessel(0))
        bfs_queue = deque()
        bfs_queue.append(head_node)
    
        # bfs traverse
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            # if it is a branch, scan until you find non-single junctions
            if cur_node.type == 'branch':
                # scan for internal junctions
                while True:
                    # retrieve outlet vessels of junc
                    next_junc = junc_mat[cur_node.ids[-1]]
                    if next_junc is None: # we reached outlet.
                        break
                    child_vess = next_junc['outlet_vessels']
                    # internal junction so add to same node
                    if len(child_vess) == 1:
                        cur_node.ids.append(child_vess[0])
                        cur_node.vessel_info.append(self.get_vessel(child_vess[0]))
                    # otherwise, exit while loop
                    else:
                        break
                if next_junc is None:
                    cur_node.children = []
                    continue
                # non-internal junction, last child, add as junction node.
                child_node = self.JunctionNode(junction_name=next_junc['junction_name'], vessel_info=next_junc)
                child_node.parent = cur_node
                cur_node.children = [child_node]
                bfs_queue.append(child_node)
            
            if cur_node.type == 'junction':
                # get outlet vessels
                child_vess = cur_node.vessel_info[-1]['outlet_vessels']
                # determine children
                cur_node.children = [self.BranchNode(vess_id=child_id, vessel_info=self.get_vessel(child_id)) for child_id in child_vess]
                for child_node in cur_node.children:
                    child_node.parent = cur_node
                    bfs_queue.append(child_node)
        
        return head_node
    
    def get_mpa_branch(self) -> BranchNode:
        ''' retrieves the MPA as a single Branch Node
        '''
        mpa = self.BranchNode(vess_id=0, 
                            vessel_info=self.get_vessel(0))
        junc_mat = self._junction_matrix()
        while True:
            child_vess = junc_mat[mpa.ids[-1]]['outlet_vessels']
            # internal junction so add to same node
            if len(child_vess) == 1:
                mpa.ids.append(child_vess[0])
                mpa.vessel_info.append(self.get_vessel(child_vess[0]))
            else:
                break
        return mpa

    def group_tree_by_generation(self, type_ = 'branch'):
        ''' groups tree by generation'''
        gens = defaultdict(list)
        
        def _dfs(root, g):
            if root.type == 'junction':
                if type_ == 'junction' or type_ == 'all':
                    gens[g].append(root)
                for child in root.children:
                    _dfs(child, g + 1)
            elif root.type == 'branch' or type_ == 'all':
                if type_ == 'branch':
                    gens[g].append(root)
                if len(root.children) == 1: # internal node, don't increase generations
                    _dfs(root.children[0], g)
            
            return

        tree = self.get_tree()
        _dfs(tree, 0)
        return gens
    
    

    def tree_bfs_iterator(self, tree, allow = "all") -> Generator[Union[BranchNode,JunctionNode], None, None]:
        ''' Tree BFS Iterator for going through the entire tree with a filter.
        '''
        if allow not in {"all", "branch", "junction"}:
            raise ValueError("Filter must be all, branch, or junction.")
        
        bfs_queue = deque()
        bfs_queue.append(tree)
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            for child_node in cur_node.children:
                bfs_queue.append(child_node)
            if allow == "all" or cur_node.type == allow:
                yield cur_node
        
