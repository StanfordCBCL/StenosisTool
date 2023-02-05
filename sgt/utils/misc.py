import os
import sys

def m2d(val, units = 'cm'):
    '''convert mmHg to dynes/cm^2 or mm^2
    '''
    if units == 'cm':
        return val * 1333.22
    elif units == 'mm':
        return val * 13.3322

def d2m(val, units = 'cm'):
    '''convert dynes/cm^2 or mm^2 to mmHg
    '''
    if units == 'cm':
        return val / 1333.22
    elif units == 'mm':
        return val / 13.3322

def blockPrint():
    ''' block printing 
    '''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    ''' enable printing
    '''
    sys.stdout = sys.__stdout__
    