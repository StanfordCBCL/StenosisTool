import os
import configparser
from collections import defaultdict
import xml.etree.ElementTree as ET
import json

class NotFound(Exception):
    pass

def check_exists(file, err = None, mkdir = False):
    ''' check if file/path exists'''
    if err == None:
        err = file + ' does not exist'
    if not os.path.exists(file):
        if mkdir:
            print(err + ': making directory')
            os.mkdir(file)
        else:
            raise NotFound(err)
    return file

def check_exists_bool(file, err = None, ignore = False):
    ''' check if file/path exists'''
    if err == None:
        err = file + ' does not exist'
    if not os.path.exists(file):
        if ignore:
            print(err)
            return False
        else:
            raise NotFound(err)
    return True


def parse_config(config_file, default_section = configparser.DEFAULTSECT):
    ''' turn a config file to a dict (all strs) '''
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
    ''' parses mdl file for faceid corresponding to caps'''
    
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

def write_json(fp, data):
    with open(fp, 'w') as sfile:
            json.dump(data, sfile, indent = 4, sort_keys=True)