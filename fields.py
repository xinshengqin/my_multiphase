#!/usr/bin/python
import numpy as np
import sys
import os
from mesh import mesh_1d,mesh_2d

class fields_2d(object):
    def __init__(self,mesh_2d, BC, IC='default', name='default_fields'):
        #BC1 is for start point
        #BC2 is for end point
        self.mesh = mesh_2d
        self.name = name
        self.BC = BC
        if IC == 'default':#initialize velocity fields as zero
            self.value=np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        else:
            self.value = IC(self.mesh)#passed from outside of class
            #user defined IC function should return a numpy array
    def apply_BC(self):
        self.value = self.BC(self.mesh,self.value)#passed from outside of class
        #para is for additional parameters

    def write_all(self,path='./'):
        file = open(path+'/'+self.name+'.csv','w')
        file.write('#No.,vedge_x,vedge_y,hedge_x,hedge_y,center_x,center_y,'+self.name+'\n')
        vx=np.ravel(self.mesh.vedge_x)
        vy=np.ravel(self.mesh.vedge_y)
        hx=np.ravel(self.mesh.hedge_x)
        hy=np.ravel(self.mesh.hedge_y)
        cx=np.ravel(self.mesh.center_x)
        cy=np.ravel(self.mesh.center_y)
        value=np.ravel(self.value)
        xy=np.vstack((vx,vy,hx,hy,cx,cy,value))
        xy=np.transpose(xy)
        np.set_printoptions(formatter={'float': lambda x: format(x, '+9.6E')})
        for i,item in enumerate(xy):
            line = '{:>9s}'.format(str(i))+','+",".join(str(item).lstrip('[').rstrip(']').split())+'\n'
            file.write(line)
        
    def write(self,path='./'):
        file = open(path+'/'+self.name+'.csv','w')
        file.write(self.name+'\n')
        value=np.ravel(self.value)
        np.set_printoptions(formatter={'float': lambda x: format(x, '+9.6E')})
        for i,item in enumerate(value):
            #line = '{:>9s}'.format(str(i))+','+",".join(str(item).lstrip('[').rstrip(']').split())+'\n'
            line = '{:>9s}'.format(str(i))+','+str(item)+'\n'
            file.write(line)
    def __str__(self):
        summary = 'mesh size:' + str( np.shape(self.mesh) ) + '\n'
        summary = summary+ 'fields size:'+str( np.shape(self.value)) +'\n'
        summary = summary+ str(self.value)

        return str(summary)

        


