import os
import configparser
from collections import defaultdict
import xml.etree.ElementTree as ET
import json
import shutil
from pathlib import Path

class NotFound(Exception):
    pass

def check_dir(dirpath: Path, mkdir = False, ignore = False):
    ''' checks if dirpath exists
    '''
    if dirpath.is_dir():
        return dirpath
    else:
        if mkdir:
            dirpath.mkdir(parents = True, exist_ok = False)
            return dirpath
        elif ignore:
            return None
        
        raise FileNotFoundError()

def check_file(filepath: Path, ignore = False):
    ''' check if file exists
    '''
    if filepath.is_file():
        return filepath
    else:
        if ignore:
            return None
        raise FileNotFoundError()

def parse_config(config_file, default_section = configparser.DEFAULTSECT):
    ''' turn a config file to a dict (all strs)
    '''
    config = configparser.ConfigParser(allow_no_value=True,
                                       default_section=default_section,
                                       interpolation=configparser.ExtendedInterpolation())
    config.read(config_file)
    
    cfg_dict = defaultdict(dict)
    for section in config.sections():
        for option in config.options(section):
            cfg_dict[section][option] = config.get(section, option)
    return cfg_dict

def parse_mdl(mdl_file, reverse = False):
    ''' parses mdl file for faceid corresponding to caps
    '''
    
    # since .mdl has a line like "<format version="1.0" />" which fails for the standard XMLparser, rather than creating a custom parser, just remove that line after reading the file in and parse as a list of strings
    with open(mdl_file, 'r') as mdl:
        lines = mdl.readlines()
        if 'format' in lines[1]:
            lines[1] = ''
    
    root = ET.fromstringlist(lines)
    faces = root[0][0].find('faces')

    # save to a dict
    face_mappings = {}
    for face in faces:
        if face.attrib['type'] == 'cap':
            if reverse:
                face_mappings[int(face.attrib['id'])] = face.attrib['name']
            else:
                face_mappings[face.attrib['name']] = int(face.attrib['id'])
    return face_mappings

def parse_face_names(datfile):
    names = []
    with open(datfile, 'r') as dfile:
        for line in dfile:
            names.append(line.rstrip())
            
    return names

def read_json(fp):
    with open(fp, 'r') as sfile:
        return json.load(sfile)
    
def write_json(fp, data):
    with open(fp, 'w') as sfile:
            json.dump(data, sfile, indent = 4, sort_keys=True)
            
def copy_rel_files(orig_dir, new_dir, exclude_solver = False):
    ''' Copies relevant files for a solver from orig dir to new dir. New dir must exist'''
    join_s = lambda file: os.path.join(orig_dir, file)
    _, dirs, files = next(os.walk(orig_dir))
    for f in files:
        fpath = join_s(f)
        if f in {'inflow.png', 'model_centerlines.vtp'}: # specific case not to be ignored
            shutil.copy(fpath,new_dir)
        elif exclude_solver and os.path.splitext(f)[-1] == '.in':
            continue
        elif os.path.splitext(f)[-1] not in {'.csv', '.npy', '.png', '.vtp', '.cfg'}:
            shutil.copy(fpath, new_dir)
        