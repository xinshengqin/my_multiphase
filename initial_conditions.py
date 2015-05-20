#!/usr/bin/python

import numpy as np

def IC_test_interpolation_u(mesh):
    value=np.array([[ 1,   2,  4,  0],
                    [ 9,   2,  3,  8],
                    [ 4,  10,  2,  2],
                    [ 3,  10,  4,  9],
                    [ 5,   9,  6,  9],
                    [ 8,   8,  4,  0]])
    return value 

def IC_test_interpolation_v(mesh):
    value=np.array([[ 3,  0,  0,  3],
                    [ 6,  4,  1,  1],
                    [ 2,  4, 10,  8],
                    [ 2, 10,  1,  8],
                    [10,  3,  1,  5],
                    [ 6, 10,  4,  9]])
    return value

def IC_test_interpolation_p(mesh):
    value=np.array([[ 2,  7,  4,  7],
                    [10,  8,  5,  2],
                    [ 2,  2,  7,  5],
                    [ 8,  3,  6,  8],
                    [ 4,  6,  6,  2],
                    [ 4,  9,  5,  7]])
    return value

def IC_flatplate_u(mesh):
    U_inf = 1.0
    value = np.zeros(mesh.center_x.shape) + U_inf
    return value

def IC_flatplate_v(mesh):
    value = np.zeros(mesh.center_x.shape) 
    return value

def IC_flatplate_p(mesh):
    value = np.zeros(mesh.center_x.shape) 
    return value
