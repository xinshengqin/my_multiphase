#!/usr/bin/python
import numpy as np
import sys
import os


def BC_fixedValue(mesh,scalar_list,value=0):
    scalar_list[:,0] = value
    scalar_list[:,-1] = value
    scalar_list[0,:] = value
    scalar_list[-1,:] = value
    return scalar_list

def BC_fixedValue_1(mesh,scalar_list):
    return BC_fixedValue(mesh,scalar_list,1)

def BC_empty(mesh,scalar_list):
    return scalar_list

def BC_zeroGradient(mesh,scalar_list):
    scalar_list[1:-1,0] = scalar_list[1:-1,1]
    scalar_list[1:-1,-1] = scalar_list[1:-1,-2]
    scalar_list[0,1:-1] = scalar_list[1,1:-1]
    scalar_list[-1,1:-1] = scalar_list[-2,1:-1]
    #need to figure out how to determine values at corner, like u[0,0]
    scalar_list[0,0] = scalar_list[0,1]
    scalar_list[0,-1] = scalar_list[0,-2]
    scalar_list[-1,0] = scalar_list[-1,1]
    scalar_list[-1,-1] = scalar_list[-1,-2]
    return scalar_list

def BC_p(mesh,scalar_list):
    
    BC_zeroGradient(mesh,scalar_list)
    scalar_list[-1,:] = 0
    return scalar_list

#BC for flow over a flat plate

def BC_flatplate_u(mesh,scalar_list):
    U_inf = 1.0
    plate_length_fraction = 0.5
    nx_before = int(round((mesh.nx+2)*plate_length_fraction))
    scalar_list[:,0] = U_inf 
    scalar_list[:,-1] = scalar_list[:,-2]
    scalar_list[0,0:nx_before] = scalar_list[1,0:nx_before]
    scalar_list[0,nx_before:] = -1*scalar_list[1,nx_before:]
    scalar_list[-1,:] = scalar_list[-2,:]
    return scalar_list

def BC_flatplate_v(mesh,scalar_list):
    plate_length_fraction = 0.5
    nx_before = int(round((mesh.nx+2)*plate_length_fraction))
    scalar_list[:,0] = 0.0
    scalar_list[:,-1] = scalar_list[:,-2]
    scalar_list[0,0:nx_before] = scalar_list[1,0:nx_before]
    scalar_list[0,nx_before:] = 0.0
    scalar_list[-1,:] = scalar_list[-2,:]
    return scalar_list
