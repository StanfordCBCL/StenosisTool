
from collections import defaultdict
import shutil
from pathlib import Path

from .misc import *
from .file_io import check_dir, parse_config

class LPNDir():
    ''' Information for an LPN Dir. Used in Tool Scripts
    '''
    
    INFLOW = 'inflow.flow'
    RCRT = 'rcrt.dat'
    CENT = 'model_centerlines.vtp'
    
    def __init__(self, lpn_dir):
        self.lpn_root = check_dir(Path(lpn_dir))
        self.lpn_path = self._find_lpn()
        self.inflow_file = self.lpn_root / self.INFLOW
        self.rcrt_file = self.lpn_root / self.RCRT
        self.centerlines = self.lpn_root / self.CENT
        self.other_base_files = {'CapInfo', 'inlet_face_names.dat', 'inlet_mapping.dat'
                                'outlet_face_names.dat', 'outlet_mapping.dat', 'rcrt.dat'}
        
        # possible dirs. They may not exist for child dirs #! correct
        self.three_d_dir = check_dir(self.lpn_root / 'three_d_dir', ignore = True)
        self.artificial_stenosis_dir = check_dir(self.lpn_root / 'artificial_sten', ignore = True)
        self.tuning_dir = check_dir(self.lpn_root / 'tuning', ignore = True)
        
    def _find_lpn(self):
        for path in self.lpn_root.iterdir():
            if path.is_file() and path.suffix == '.in':
                return path
        raise FileNotFoundError('LPN file not found. Please set the suffix to .in')

    

class ModelPath():
    ''' Model specific info & Paths
    '''
    
    CONFIG_DIR = 'config_files'
    CENTERLINES_DIR = 'centerline_files'
    BASE_LPN_DIR = 'base_lpn_files'
    TUNING_DIR = 'tuning_files'
    
    MODEL_CENT = 'model_centerlines.vtp'
    
    MODEL_LPN = 'lpn.in'
    TUNING_LPN = 'tuning.in'
    
    STENOSIS_FILE = 'stenosis_file.txt'
    DEV_CONF = 'dev_config.ini'


    def __init__(self, root : Path, model_type: str):
        self.model_root = root
        self.type = model_type
        self.model_name = self.model_root.name
        
        # retrieve dev config
        self.config_files = check_dir(self.model_root / self.CONFIG_DIR, mkdir = True)
        
        # search and read for dev_config file in self.config_files
        self.dev_config = self.config_files / self.DEV_CONF
        if not self.dev_config.is_file():
            self.construct_dev_config(self.dev_config)
            
        # read dev config
        self.info = parse_config(str(self.dev_config))
        
        # check info for completeness
        if not self.check_info():
            print('Info for model {} is unavailable'.format(self.model_name))
            return
        
        # joins root to the files as a Path object.
        for key, vals in self.info['files'].items():
            if vals != '':
                self.info['files'][key] = str(self.model_root / vals)
        
        # dirs            
        self.centerline_dir = check_dir(self.model_root / self.CENTERLINES_DIR, mkdir = True)
        self.base_lpn_dir = check_dir(self.model_root / self.BASE_LPN_DIR, mkdir = True)
        
        # files
        self.base_lpn = self.base_lpn_dir / (self.info['metadata']['name'] + '_' + self.MODEL_LPN)
        self.model_centerlines = self.centerline_dir / (self.info['metadata']['name'] + '_' + self.MODEL_CENT)

    def check_info(self):
        ''' validates metadata exists
        '''
        for data in self.info['metadata']:
            if data == '':
                return False
        return True

    def construct_dev_config(self, config_fp: Path):
        ''' constructs an empty dev config file from template
        '''
        template_config = Path(__file__).parent / 'templates/dev_config.ini'
        shutil.copy(str(template_config), str(config_fp))
        
    def __repr__(self):
        s = ''
        for attr, val in self.__dict__.items():
            s += '{} : \n\t{}\n'.format(attr, val)
        return s
        
class DataPath():
    ''' Data Manager
    '''
    
    DATA = 'data'
    HEALTHY = 'healthy'
    STENOSIS = 'stenosis'
        
    def __init__(self, root):
        # main data dirs
        self.project_root = check_dir(Path(root))
        self.data_root = check_dir(self.project_root / self.DATA, mkdir = True)
        self.healthy = check_dir(self.data_root / self.HEALTHY, mkdir = True)
        self.stenosis = check_dir(self.data_root / self.STENOSIS, mkdir = True)
        
        # retrieve models
        self.models = defaultdict(dict)
        self.model_names = set()
        for model_path in self.healthy.iterdir():
            if not model_path.is_dir():
                pass
            self.models[self.HEALTHY][model_path.name] = ModelPath(model_path, self.HEALTHY)
            self.model_names.add(model_path.name)
                
        for model_path in self.stenosis.iterdir():
            if not model_path.is_dir():
                pass
            self.models[self.STENOSIS][model_path.name] = ModelPath(model_path, self.STENOSIS)
            self.model_names.add(model_path.name)
                
    def __repr__(self):
        s = ''
        for attr, val in self.__dict__.items():
            s += '{} : \n\t{}\n'.format(attr, val)
        return s[1:]

    def run_all(self, func, *args):
        ''' wrapper for functions utilizing all models. Function must reserve first position for a ModelPath object
        '''
        output = {}
        for type, modeldict in self.models.items():
            for mod_name, mod_path in modeldict.items():
                print('Running function for {}...'.format(mod_name))
                output[mod_name] = func(mod_path, *args)
        return output

    def find_model(self, model_name) -> ModelPath:
        ''' finds a model
        '''
        if model_name in self.models[self.HEALTHY]:
            return self.models[self.HEALTHY][model_name]
        elif model_name in self.models[self.STENOSIS]:
            return self.models[self.STENOSIS][model_name]
        else:
            print(model_name + ' is an invalid model')
            return None

