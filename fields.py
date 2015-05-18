#!/usr/bin/python
import numpy as np
import sys
import os
from mesh import Mesh_1d,Mesh_2d

class Fields_2d(object):
    def __init__(self,mesh_2d, BC =None, IC='default',ftype='center', name='default_fields'):
        self.mesh = mesh_2d
        self.name = name
        self.BC = BC
        self.IC = IC
        self.ftype = ftype
        self.center= None
        self.corner= None
        self.vedge= None
        self.hedge= None
        self.fdict = {'center':self.center,'corner':self.corner, 'vedge':self.vedge,'hedge':self.hedge}
        if IC == 'default':#initialize velocity fields as zero
            self.fdict[ftype] = np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        else:
            #todo IC should take more arguments
            self.fdict[ftype] = IC(self.mesh)#passed from outside of class
            #user defined IC function should return a numpy array
        self.value = self.fdict[self.ftype]
    def applyBC(self):
        self.fdict[self.ftype] = self.BC(self.mesh,self.value)#passed from outside of class
        self.value = self.fdict[self.ftype]
        #para is for additional parameters

    def compute_center(self):
        self.center = np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        if self.ftype == 'center':
            self.center = self.value
        elif self.ftype == 'hedge':
            for j in range(1,self.mesh.ny+4):
                self.center[j,:] = (self.value[j,:]+self.value[j-1,:])/2.
        elif self.ftype == 'vedge':
            for i in range(1,self.mesh.nx+4):
                self.center[:,i] = (self.value[:,i]+self.value[:,i-1])/2.
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        return self.center

    def compute_vedge(self):
        self.vedge = np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        if self.ftype == 'center':
            for i in range(0,self.mesh.nx+3):
                self.vedge[:,i] = (self.value[:,i+1]+self.value[:,i])/2.
        elif self.ftype == 'hedge':
            for i in range(0,self.mesh.nx+3):
                for j in range(1,self.mesh.ny+4):
                    self.vedge[j,i] = (self.value[j-1,i]+self.value[j-1,i+1]+self.value[j,i]+self.value[j,i+1])/4.
        elif self.ftype == 'vedge':
            self.vedge = self.value
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        return self.vedge

    def compute_hedge(self):
        self.hedge = np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        if self.ftype == 'center':
            for j in range(0,self.mesh.ny+3):
                self.hedge[j,:] = (self.value[j,:]+self.value[j+1,:])/2.
        elif self.ftype == 'hedge':
            self.hedge = self.value
        elif self.ftype == 'vedge':
            for i in range(1,self.mesh.nx+4):
                for j in range(0,self.mesh.ny+3):
                    self.hedge[j,i] = (self.value[j,i]+self.value[j,i-1]+self.value[j+1,i]+self.value[j+1,i-1])/4.
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        return self.hedge

    def compute_corner(self):
        self.corner = np.zeros((self.mesh.ny+4,self.mesh.nx+4))
        if self.ftype == 'center':
            for j in range(0,self.mesh.ny+3):
                for i in range(0,self.mesh.nx+3):
                    self.corner[j,i] = (self.value[j,i]+self.value[j,i+1]+self.value[j+1,i]+self.value[j+1,i+1])/4.
        elif self.ftype == 'hedge':
            for i in range(0,self.mesh.nx+3):
                self.corner[:,i] = (self.value[:,i+1]+self.value[:,i])/2.
        elif self.ftype == 'vedge':
            for j in range(0,self.mesh.ny+3):
                self.corner[j,:] = (self.value[j+1,:]+self.value[j,:])/2.
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        return self.corner

    def ddx2(self):#d2(self)/dx2
        ddx2 = Fields_2d(self.mesh,ftype = self.ftype)
        for i in range(1,self.mesh.nx+3):
            ddx2.value[:,i] = (self.value[:,i+1]-2*self.value[:,i]+self.value[:,i-1])/(self.mesh.dx**2.)
        return ddx2
    def ddy2(self):#d2(self)/dy2
        ddy2 = Fields_2d(self.mesh,ftype = self.ftype)
        for j in range(1,self.mesh.ny+3):
            ddy2.value[j,:] = (self.value[j+1,:]-2*self.value[j,:]+self.value[j-1,:])/(self.mesh.dy**2.)
        return ddy2
    
    def __add__(self,fields):
        #todo add a function to check if two fields can be added
        if (self.ftype != fields.ftype): 
            print "try to add up:",self.name,"and",fields.name
            print "cannot add two fields of different types up!"
            sys.exit()
        else:
            result = Fields_2d(self.mesh,ftype=self.ftype)
            result.value = self.value + fields.value
            return result

    def __mul__(self,scalar):
        if type(scalar) != float:
            print "multiply operation must be taken between floats and fields"
            sys.exit()
        else:
            result = Fields_2d(self.mesh,ftype=self.ftype)
            result.value = self.value * scalar
            return result
    def __rmul__(self,scalar):
        if type(scalar) != float:
            print "multiply operation must be taken between floats and fields"
            sys.exit()
        else:
            result = Fields_2d(self.mesh,ftype=self.ftype)
            result.value = self.value * scalar
            return result

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
        if path != './':
            if not( path in os.listdir('./')):
                os.mkdir( path )
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

        


