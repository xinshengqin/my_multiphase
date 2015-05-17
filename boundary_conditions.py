#!/usr/bin/python
import numpy as np
import sys
import os


def BC_fixedValue(mesh,scalar_list,value=0):
    scalar_list[:,0] = value
    scalar_list[:,1] = value
    scalar_list[:,-1] = value
    scalar_list[:,-2] = value
    scalar_list[0,:] = value
    scalar_list[1,:] = value
    scalar_list[-1,:] = value
    scalar_list[-2,:] = value
    return scalar_list

def BC_fixedValue_1(mesh,scalar_list):
    return BC_fixedValue(mesh,scalar_list,1)
