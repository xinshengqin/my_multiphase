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

def IC_C_circle(mesh):
    #initial condition for volume fraction C
    #C is 1 in a circle centered at (5,5) with a radius of 4
    value = np.zeros(mesh.center_x.shape) 
    resolution = 20
    xc = 5.
    yc = 5.
    radius = 4
    for j in range(1,mesh.ny+1):
        for i in range(1,mesh.nx+1):
            area = 0.
            x = np.linspace(mesh.vedge_x[j,i-1],mesh.vedge_x[j,i],resolution)
            y = np.linspace(mesh.hedge_y[j-1,i],mesh.hedge_y[j,i],resolution)
            x = x[1:-1]
            y = y[1:-1]
            #compute C in cell[j,i]
            for k in range(x.size):
                for l in range(y.size):
                    if ( (x[k]-xc)**2+(y[l]-yc)**2 ) < radius**2:
                        area = area+mesh.dx*mesh.dy/(x.size*y.size)
            value[j,i] = area/(mesh.dx*mesh.dy)

    return value



