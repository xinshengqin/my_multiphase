#!/usr/bin/python
import numpy as np
import sys
import os
from mesh import mesh_1d,mesh_2d

class fields_2d(object):
    def __init__(self,mesh_2d,IC='default', BC = 'periodic', name='default_fields'):
        self.mesh = mesh_2d
        self.name = name
        if IC == 'default':#initialize velocity fields as zero
            self.scalar=np.zeros(self.mesh.nx+4,self.mesh.ny+4)
            self.v=np.zeros(self.mesh.nx+4,self.mesh.ny+4)
            self.p=np.zeros(self.mesh.nx+4,self.mesh.ny+4)
        else:
            [self.u,self.v,self.p] = IC(self.mesh)#passed from outside of class
            #user defined IC function should return u,v,p

        #self.BC=BC
        #if BC['type'] == 'dirichlet':
        #    pass#todo
        #elif BC['type'] == 'neumann':
        #    pass#todo
        #elif BC['type'] == 'periodic':
        #    self.BC_periodic(self.u)
        #    self.BC_periodic(self.v)
        #    self.BC_periodic(self.p)
        #else:
        #    print 'BC type cannot found in database'
        #    sys.exit()
    def apply_BC(self,BC):
        [self.u,self.v,self.p] = BC(self.mesh,self.u,self.v,self.p)#passed from outside of class
        

    def BC_periodic(self,scalar_list):
        scalar_list[2] = scalar_list[-2]
        scalar_list[3] = scalar_list[-1]
        scalar_list[-4] = scalar_list[0]
        scalar_list[-3] = scalar_list[1]

