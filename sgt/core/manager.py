# File: manager.py
# File Created: Monday, 31st October 2022 6:02:27 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Saturday, 3rd December 2022 1:14:51 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Manages files for 1 model in pipeline



from pathlib import Path
import shutil

from sgt.utils import io

class Manager():
    ''' Base Manager Class
    '''
    
    def __init__(self, config: str):
        self.config = Path(config)
        files = io.read_json(self.config)
        self.info = files
        
        # Options
        self.tune = files['options']['tune']
        self.jc = files['options']['jc'] if 'jc' in files['options'] else False
        self.lpn_suff = '.jc.in' if self.jc else '.in'
        
        # File Paths
        self.root = io.check_dir(Path(files['workspace']['root']))
        self.surface_model = io.check_file(self.root / files['workspace']['surface_model'])
        self.centerlines = io.check_file(self.root / files['workspace']['centerlines'])
        self.mdl = io.check_file(self.root / files['workspace']['mdl_mapping'])
        self.flow = io.check_file(self.root / files['workspace']['flow_file'])
        self.rcrt = io.check_file(self.root / files['workspace']['rcrt_file'], ignore = True)
        
        self.capinfo = io.check_file(self.root / files['workspace']['capinfo']) if self.tune else io.check_file(self.root / files['workspace']['capinfo'], ignore = True)
        # MetaData
        self.model_name = files['metadata']['model_name']
        self.diseased = files['metadata']['diseased']
        self.inlet = files['metadata']['inlet']
        self.units = files['metadata']['units'] if 'units' in files['metadata'] else 'cm'
        
        # tune params
        self.tune_params = files['tune_params'] if self.tune else None
        if 'tune_params' not in files and self.tune:
            raise ValueError("tune_params not found despite tune = True")
        
        
        
        
    
    @staticmethod
    def construct_dev_config(config_fp: Path):
        ''' constructs an empty config json from template
        '''
        template_config = Path(__file__).parent.parent / 'templates/config.json'
        if (config_fp / 'config.json').exists():
            print("Config File already exists.")
            return
        shutil.copy(str(template_config), str(config_fp))
        print("Config Written.")
        
    
    def config_add(self, indices: list, value):
        ''' add a value to config
        '''
        cur_layer = self.info
        for i in indices[:-1]:
            if i not in cur_layer:
                cur_layer[i] = {}
            cur_layer = cur_layer[i]
        cur_layer[indices[-1]] = value
    
    def write_config(self):
        ''' writes to config
        '''
        io.write_json(self.config, self.info)
            
    
class SVManager(Manager):
    ''' Manages file in SV directory form 
    '''

    def __init__(self, config):
        super().__init__(config)
    
    def construct_workspace(self, outdir: Path, force = True):
        ''' constructs a workspace for the pipeline by copying files over
        '''
        if outdir.exists() and not force:
            raise FileExistsError(outdir)
        elif outdir.exists():
            shutil.rmtree(str(outdir))
            
        outdir.mkdir(parents = True, exist_ok = True)
        
        ## Copy
        # surface model
        sm = outdir / (self.model_name + '.vtp')
        shutil.copy(str(self.surface_model), str(sm))
        
        # centerlines
        c = outdir / (self.model_name + '_centerlines.vtp')
        shutil.copy(str(self.centerlines), str(c))
        
        # mdl
        mdl = outdir / (self.model_name + '.mdl')
        shutil.copy(str(self.mdl), str(mdl))
        
        # flow
        flow = outdir / "inflow.flow"
        shutil.copy(str(self.flow), str(flow))
        
        # BC if precomputed
        if self.rcrt:
            bc = outdir / self.rcrt.name
            shutil.copy(str(self.rcrt), str(bc))
            
        if self.tune:
            capinfo = outdir / "CapInfo"
            shutil.copy(str(self.capinfo), str(capinfo) )
        
        # write a config json
        cfg = {
            'workspace': {
                'root' : str(outdir),
                'surface_model' : self.model_name + '.vtp',
                'mdl_mapping' : self.model_name + '.mdl',
                'centerlines' : self.model_name + '_centerlines.vtp',
                'flow_file' :  "inflow.flow",
                'rcrt_file' : str(self.rcrt.name) if self.rcrt else '',
                'capinfo' : "CapInfo" if self.capinfo else ''
            },
            'metadata': {
                'model_name' : self.model_name,
                'diseased' : self.diseased,
                'inlet' : self.inlet
            },
            'options': {
                'tune': self.tune
            }
        }
        
        # add rest of meta data in.
        for key, val in self.info['metadata'].items():
            if key not in cfg['metadata']:
                cfg['metadata'][key] = val
        
        # if tune
        if self.tune:
            cfg['tune_params'] = self.tune_params
        
        config = outdir / 'config.json'
        io.write_json(config, cfg)
        
class LPNConstructionManager(Manager):
    ''' manages LPN creation 
    '''
    def __init__(self, config: str):
        super().__init__(config)
        
        
        self.lpn_files = io.check_dir(self.root / "LPN_files", mkdir = True)
        self.lpn = self.lpn_files / (self.model_name + self.lpn_suff)
        
        

class TuningManager(LPNConstructionManager):
    ''' manages Tuning
    '''
    
    def __init__(self, config: str):
        super().__init__(config)
        
        # tuning dir
        self.tuning_dir = io.check_dir(self.lpn_files / "tuning", mkdir = True)
        
        # files
        self.tuning_lpn = self.tuning_dir / (self.model_name + '_tuning.in')
        self.results = self.tuning_dir / ('results.npy')
        
        # intermediate dirs
        self.intermediate_dir = io.check_dir(self.tuning_dir / "intermediate_tuning", mkdir = True)
        
        # sensitivity test dir
        self.sens_dir = io.check_dir(self.tuning_dir / "sens_test", mkdir = True)
        
class StenosisParametrizationManager(LPNConstructionManager):
    ''' manages artificial stenosis
    '''
    
    def __init__(self, config):
        super().__init__(config)
        
        if self.diseased: # stenosis
            
            # real stenosis
            self.SP_dir = io.check_dir(self.lpn_files / 'real_stenosis', mkdir = True)
            
            # copy of fixed stenosis lpn
            self.SP_lpn = self.SP_dir / (self.lpn.stem.split('.')[0] + '.fix' + ''.join(self.lpn.suffixes))
            
        else: # healthy
            
            # artificial stenosis dir.
            self.SP_dir = io.check_dir(self.lpn_files / 'artificial_stenosis', mkdir = True)
            
            # new artificial stenosis lpn
            self.SP_lpn = self.SP_dir / (self.lpn.stem.split('.')[0] + '.as' + ''.join(self.lpn.suffixes))
            
            # if a config file to output
            self.SP_config = self.SP_dir / 'AS_config.json'
            
        # changes
        self.stenosis_parametrization = self.SP_dir / 'stenosis_parametrization.dat'
        

class DataGenerationManager(StenosisParametrizationManager):
    
    def __init__(self, config):
        super().__init__(config)
        
        self.data_dir = io.check_dir(self.SP_dir / 'data', mkdir = True)
        
        # select the appropriate lpn to use as the data generation base.
        if self.diseased:
            self.data_gen_lpn = self.lpn
        else:
            self.data_gen_lpn = self.SP_lpn