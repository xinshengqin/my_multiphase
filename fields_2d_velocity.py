#!/usr/bin/python
import numpy as np
import sys
import os
from fields import Fields_2d

class U(Fields_2d):
    def __init__(self,mesh_2d, BC = None, IC = 'default', name = 'U'):
        Fields_2d.__init__(self,mesh_2d,BC,IC, ftype = 'vedge',name = 'U')
    def uphi_x_vedge(self,phi):#\frac{\partial u*phi}{\partial x} on vertical edge
        #phi is prefered to be defined on vertical edge
        uphix = Fields_2d(self.mesh,ftype = 'vedge')
        self.compute_center()
        phi.compute_center()
        for i in range(0,self.mesh.nx+3):
            uphix.value[:,i] = 1./self.mesh.dx * (self.center[:,i+1]*phi.center[:,i+1]-self.center[:,i]*phi.center[:,i])
        return uphix

    def uphi_x_hedge(self,phi):#\frac{\partial u*phi}{\partial x} on horizontal edge
        #phi prefered to be defined on horizontal edge
        uphix = Fields_2d(self.mesh,ftype = 'hedge')
        self.compute_corner()
        phi.compute_corner()
        for i in range(1,self.mesh.nx+4):
            uphix.value[:,i] = 1./self.mesh.dx * (self.corner[:,i]*phi.corner[:,i]-self.corner[:,i-1]*phi.corner[:,i-1])
        return uphix






class V(Fields_2d):
    def __init__(self,mesh_2d, BC = None, IC = 'default', name = 'V'):
        Fields_2d.__init__(self,mesh_2d,BC, IC, ftype = 'hedge',name = 'V')
        self.BC = BC

    def vphi_y_vedge(self,phi):#\frac{\partial v*phi}{\partial y} on vertical edge
        vphiy = Fields_2d(self.mesh,ftype='vedge')
        self.compute_corner()
        phi.compute_corner()
        for j in range(1,self.mesh.ny+4):
            vphiy.value[j,:] = 1./self.mesh.dy * (self.corner[j,:]*phi.corner[j,:]-self.corner[j-1,:]*phi.corner[j-1,:])
        return vphiy


    def vphi_y_hedge(self,phi):#\frac{\partial v*phi}{\partial y} on vertical edge
        vphiy = Fields_2d(self.mesh,ftype='hedge')
        self.compute_center()
        phi.compute_center()
        for j in range(0,self.mesh.ny+3):
            vphiy.value[j,:] = 1./self.mesh.dy * (self.center[j+1,:]*phi.center[j+1,:]-self.center[j,:]*phi.center[j,:])
        return vphiy

