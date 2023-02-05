# File: io.py
# File Created: Monday, 31st October 2022 6:06:44 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 3rd November 2022 9:29:03 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: IO utilities
import json
from pathlib import Path
import xml.etree.ElementTree as ET

def read_json(fp: Path):
    ''' reads a json as dict
    '''
    fp = str(fp)
    with open(fp, 'r') as sfile:
        return json.load(sfile)
    
def write_json(fp: Path, data, indent = 4, sort_keys = True):
    ''' writes a dict as json
    '''
    fp = str(fp)
    with open(fp, 'w') as sfile:
            json.dump(data, sfile, indent = indent, sort_keys=sort_keys)


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
        
        raise FileNotFoundError(dirpath)

def check_file(filepath: Path, ignore = False):
    ''' check if file exists
    '''
    if filepath.is_file():
        return filepath
    else:
        if ignore:
            return None
        raise FileNotFoundError(filepath)
    
def parse_mdl(mdl_file: Path, reverse = False):
    ''' parses mdl file for faceid corresponding to caps
    '''
    
    # since .mdl has a line like "<format version="1.0" />" which fails for the standard XMLparser, rather than creating a custom parser, just remove that line after reading the file in and parse as a list of strings
    mdl_file = str(mdl_file)
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
