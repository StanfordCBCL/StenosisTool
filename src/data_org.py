
from collections import defaultdict
from os import path
import os

from .misc import *
from .file_io import check_exists, check_exists_bool, parse_config

class DataPath():
    ''' Data Manager '''
    
    HEALTHY = 'healthy'
    STENOSIS = 'stenosis'
    
    class ModelPath():
        ''' Model specific info '''
        
        CONFIG_DIR = 'config_files'
        CENTERLINES_DIR = 'centerline_files'
        SOLVER_DIR = 'solver_files'
        ARTIFICIAL_STEN_DIR = 'artificial_stenosis'
        FIXED_STEN_DIR = 'fixed_stenosis_dir'
        
        MODEL_CENT = 'model_centerlines.vtp'
        TUNING_CENT = 'tuning_centerlines.vtp'
        
        MODEL_SOLVER = 'model.in'
        TUNING_SOLVER = 'tuning.in'
        
        STENOSIS_FILE = 'stenosis_file.txt'
        DEV_CONF = 'dev_config.ini'
    
    
        def __init__(self, root, type):
            self.model_root = root
            
            # config
            self.config_files, _ = check_exists(path.join(root, self.CONFIG_DIR), 
                                                err = 'Model {} does not have a {} directory'.format(path.basename(root), self.CONFIG_DIR), mkdir = True)
            
            # search and read for dev_config file in self.config_files
            self.dev_config = path.join(self.config_files, self.DEV_CONF)
            exists = check_exists_bool(self.dev_config,
                                err = 'Model {} does not have a {} file. Creating one.'.format(os.path.basename(root), self.DEV_CONF), ignore=True)
            if not exists:
                self.construct_dev_config()
                
            # read dev config
            self.info = parse_config(self.dev_config)
            for key, vals in self.params['files'].items():
                if vals != '':
                    self.params['files'][key] = path.join(root, vals)
            
            
            # dirs            
            self.centerline_dir = check_exists(path.join(root, self.CENTERLINES_DIR), mkdir = True)
            self.solver_dir = check_exists(path.join(root, self.SOLVER_DIR), mkdir = True)
            
            # files
            self.model_solver = path.join(self.solver_dir, self.info['metadata']['name'] + '_' + self.MODEL_SOLVER)
            self.tune_solver = path.join(self.solver_dir, self.info['metadata']['name'] + '_' + self.TUNING_SOLVER)
            self.model_centerlines = path.join(self.centerline_dir, self.info['metadata']['name'] + '_' + self.MODEL_CENT)
            self.tune_centerlines = path.join(self.centerline_dir, self.info['metadata']['name'] + '_' + self.TUNING_CENT)
            
            # model specific dirs
            self.type = type
            self.artificial_sten_dir = None
            self.fixed_sten_dir = None
            self.stenosis_file = None
            if type == 'healthy':
                self.artificial_sten_dir = check_exists(path.join(root, self.ARTIFICIAL_STEN_DIR), mkdir = True)
            elif type == 'stenosis':
                self.stenosis_vessels_file = path.join(self.solver_dir, self.STENOSIS_FILE)
                self.stenosis_check_dir = check_exists(path.join(root, self.FIXED_STEN_DIR), mkdir = True)
                
        def construct_dev_config(self):
            ''' constructs an empty dev config '''
            pass
        
        def __repr__(self):
            s = ''
            for attr, val in self.__dict__.items():
                s += '{} : \n\t{}\n'.format(attr, val)
            return s
        
    def __init__(self, root):
        self.data_root = root
        self.healthy = check_exists(path.join(root, self.HEALTHY), mkdir = True)
        self.stenosis = check_exists(path.join(root, self.STENOSIS), mkdir = True)
        self.data_dirs = [self.healthy, self.stenosis]
        
        # retrieve models
        self.models = defaultdict(dict)
        for mod in os.listdir(self.healthy):
            model_path = path.join(self.healthy, mod)
            if os.path.isdir(model_path):
                self.models[self.HEALTHY][mod] = self.ModelPath(model_path, type = self.HEALTHY)
                
        for mod in os.listdir(self.stenosis):
            model_path = path.join(self.stenosis, mod)
            if os.path.isdir(model_path):
                self.models[self.STENOSIS][mod] = self.ModelPath(model_path, type = self.STENOSIS)
                
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
    