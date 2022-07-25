import os
import argparse
import sys

def m2d(val):
    '''convert mmHg to dynes/cm^2'''
    return val * 1333.22

def d2m(val):
    '''convert dynes/cm^2 to mmHg'''
    return val / 1333.22

def blockPrint():
    ''' block printing '''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    ''' enable printing'''
    sys.stdout = sys.__stdout__
    
def create_parser(desc = '') -> argparse.ArgumentParser:
    ''' creates a parse for scripts '''
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(help='Mode for script', dest = 'mode')
    tool = subparsers.add_parser(name='tool')
    dev = subparsers.add_parser(name='dev')
    dev.add_argument('-root', default='.', help='Root of entire project')
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    


    return parser, dev, tool

def get_basename(fp):
    return os.path.splitext(os.path.basename(fp))[0]

# Not in use
'''
def parse_bsa_chart(bsa_file):
    with open(bsa_file, 'r') as bsafile:
        bsafile.readline()
        
        bsa_chart = {}
        for line in bsafile:
            line = line.rstrip().split()
            age = line[0]
            bsa = float(line[1])
            
            if '-' in age:
                age = age.split('-')
                for a in range(int(age[0]), int(age[1]) + 1):
                    bsa_chart[a] = bsa
            elif '+' in age:
                age = [age.split('+')[0] , 120] # 120 is pretty old
                for a in range(int(age[0]), int(age[1]) + 1):
                    bsa_chart[a] = bsa
            else:
                bsa_chart[int(age)] = bsa
                
    return bsa_chart


def parse_height_chart(height_file):
    height_chart = {}
    with open(height_file, 'r') as hfile:
        headers = hfile.readline().split('\t') # ignore
        for line in hfile:
            line = line.rstrip().split('\t')
            age = int(float(line[0]))
            median_height = line[6]
            if age / 12  == age // 12:
                height_chart[int(age//12)] = median_height
            if age / 12 == 20:
                for i in range(21, 120): # fill in up to 120
                    height_chart[i] = median_height
                    
    return height_chart
'''