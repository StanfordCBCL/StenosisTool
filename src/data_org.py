
from collections import defaultdict
from os import path
import os
import shutil

from .misc import *
from .file_io import check_exists, check_exists_bool, parse_config

class ModelPath():
    ''' Model specific info '''
    
    CONFIG_DIR = 'config_files'
    CENTERLINES_DIR = 'centerline_files'
    SOLVER_DIR = 'solver_files'
    TUNING_DIR = 'tuning_files'
    
    MODEL_CENT = 'model_centerlines.vtp'
    
    MODEL_SOLVER = 'model.in'
    TUNING_SOLVER = 'tuning.in'
    
    STENOSIS_FILE = 'stenosis_file.txt'
    DEV_CONF = 'dev_config.ini'


    def __init__(self, root, model_type):
        self.model_root = root
        self.type = model_type
        self.model_name = path.basename(root)
        
        # config
        self.config_files = check_exists(path.join(root, self.CONFIG_DIR), 
                                            err = 'Model {} does not have a {} directory'.format(path.basename(root), self.CONFIG_DIR), mkdir = True)
        
        # search and read for dev_config file in self.config_files
        self.dev_config = path.join(self.config_files, self.DEV_CONF)
        exists = check_exists_bool(self.dev_config,
                            err = 'Model {} does not have a {} file. Creating one.'.format(os.path.basename(root), self.DEV_CONF), ignore=True)
        if not exists: # exit after constructing dev config
            self.construct_dev_config(self.dev_config)
            

            
            
        # read dev config
        self.info = parse_config(self.dev_config)
        
        if not self.check_info():
            print('Info for model {} is unavailable'.format(self.model_name))
            return
        
        for key, vals in self.info['files'].items():
            if vals != '':
                self.info['files'][key] = path.join(root, vals)
        
        
        # dirs            
        self.centerline_dir = check_exists(path.join(root, self.CENTERLINES_DIR), mkdir = True)
        self.solver_dir = check_exists(path.join(root, self.SOLVER_DIR), mkdir = True)
        self.tuning_dir = check_exists(path.join(root, self.TUNING_DIR), mkdir = True)
        
        # files
        self.model_solver = path.join(self.solver_dir, self.info['metadata']['name'] + '_' + self.MODEL_SOLVER)
        self.tune_solver = path.join(self.tuning_dir, self.info['metadata']['name'] + '_' + self.TUNING_SOLVER)
        self.model_centerlines = path.join(self.centerline_dir, self.info['metadata']['name'] + '_' + self.MODEL_CENT)
    
    def check_info(self):
        for data in self.info['metadata']:
            if data == '':
                return False
        return True

    def construct_dev_config(self, dev_config_fp):
        ''' constructs an empty dev config '''
        file = '# general model information\n'
        file += '[metadata]\nid = \nname = \nage = \ngender = \ncondition = \n\n'
        
        file += '# file information\n'
        file += '[files]\n# flow files\ninflow = \nrom_inflow = \n\n# model files\nmdl_file = \nvtp_file = \ncap_info = \n'
        
        file += '# 3D model info\n'
        file += '[model]\ninlet = \nlpa = \nrpa = \nunits = \n'
    
        with open(dev_config_fp, 'w') as dfile:
            dfile.write(file)
        
    
    def __repr__(self):
        s = ''
        for attr, val in self.__dict__.items():
            s += '{} : \n\t{}\n'.format(attr, val)
        return s
        
class DataPath():
    ''' Data Manager '''
    
    HEALTHY = 'healthy'
    STENOSIS = 'stenosis'
        
    def __init__(self, root):
        self.project_root = root
        self.data_root = check_exists(path.join(root, 'data'), err= 'Data Dir does not exist')
        self.healthy = check_exists(path.join(self.data_root, self.HEALTHY), mkdir = True)
        self.stenosis = check_exists(path.join(self.data_root, self.STENOSIS), mkdir = True)
        self.data_dirs = [self.healthy, self.stenosis]
        
        # retrieve models
        self.models = defaultdict(dict)
        self.model_names = set()
        for mod in os.listdir(self.healthy):
            model_path = path.join(self.healthy, mod)
            if os.path.isdir(model_path):
                self.models[self.HEALTHY][mod] = ModelPath(model_path, self.HEALTHY)
                self.model_names.add(mod)
                
        for mod in os.listdir(self.stenosis):
            model_path = path.join(self.stenosis, mod)
            if os.path.isdir(model_path):
                self.models[self.STENOSIS][mod] = ModelPath(model_path, self.STENOSIS)
                self.model_names.add(mod)
                
    def __repr__(self):
        s = ''
        for attr, val in self.__dict__.items():
            s += '{} : \n\t{}\n'.format(attr, val)
        return s[1:]


    def run_all(self, func, *args):
        ''' wrapper for functions utilizing all models. Function must reserve first position for a ModelPath object'''
        output = {}
        for type, modeldict in self.models.items():
            for mod_name, mod_path in modeldict.items():
                print('Running function for {}...'.format(mod_name))
                output[mod_name] = func(mod_path, *args)
        return output

    def find_model(self, model_name) -> ModelPath:
        ''' finds a model'''
        if model_name in self.models[self.HEALTHY]:
            return self.models[self.HEALTHY][model_name]
        elif model_name in self.models[self.STENOSIS]:
            return self.models[self.STENOSIS][model_name]
        else:
            print(model_name + ' is an invalid model')
            return None
    
    
class Results():
    
    RESULTS_DIR = 'results'
    
    def __init__(self, root, model: ModelPath):
        
        # convert Model path to an equivalent results path
        self.root = root
        
        self.model = model
        self.results_dir = check_exists(os.path.join(self.root, self.RESULTS_DIR), mkdir = True)
        self.type_dir = check_exists(os.path.join(self.results_dir, self.model.type), mkdir = True)
        

        self.model_dir  = check_exists(os.path.join(self.results_dir, model.info['metadata']['name']), mkdir = True)

        

class StenosisToolResults(Results):
    
    BASE_SOLVER_DIR = 'base_solver_dir'
    ART_STEN_DIR = 'artificial_stenosis_dir'
    FIXED_STEN_DIR = 'fixed_stenosis_dir'
    
    def __init__(self, root, model: ModelPath):
        super().__init__(root, model)
        self.base_solver_dir = check_exists(os.path.join(self.model_dir, self.BASE_SOLVER_DIR), mkdir = True)

        self.base_solver = os.path.join(self.base_solver_dir, os.path.basename(model.model_solver))
        if not os.path.exists(self.base_solver):
            shutil.copy(model.model_solver, self.base_solver)
        
        if self.model.type == 'healthy':
            self.setup_healthy()
        elif self.model.type == 'stenosis':
            self.setup_stenosis()
            
        
        
    def setup_stenosis(self):
        ''' if self.type = None, set up like a stenosis model '''
        self.fixed_stenosis = check_exists(os.path.join(self.model_dir, self.FIXED_STEN_DIR), mkdir = True)
    
    def setup_healthy(self):
        ''' if self.type = healthy'''
        self.artificial_sten_dir = check_exists(os.path.join(self.model_dir, self.ART_STEN_DIR), mkdir = True)
        
        
        
    