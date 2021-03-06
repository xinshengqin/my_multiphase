#!/usr/bin/python
import numpy as np
import sys
import os
import copy

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
        self.fdict[ftype] = np.zeros(self.mesh.center_x.shape)

        #initialize
        if IC != 'default':#initialize velocity fields as zero
            self.fdict[ftype] = IC(self.mesh)#passed from outside of class
            #user defined IC function should return a numpy array
        self.value = self.fdict[self.ftype]
        if self.value.shape != self.mesh.center_x.shape:
            print "ERROR! The IC function return a field that does not match the size of mesh!"
            sys.exit()
    def applyBC(self):
        self.fdict[self.ftype] = self.BC(self.mesh,self.value)#passed from outside of class
        self.value = self.fdict[self.ftype]
        #para is for additional parameters

    def compute_center(self):
        center = Fields_2d(self.mesh,ftype = 'center')
        #todo may not need below in the future
        self.center = np.zeros(self.mesh.center_x.shape)
        if self.ftype == 'center':
            self.center = self.value
        elif self.ftype == 'hedge':
            for j in range(1,self.mesh.ny+2):
                self.center[j,:] = (self.value[j,:]+self.value[j-1,:])/2.
        elif self.ftype == 'vedge':
            for i in range(1,self.mesh.nx+2):
                self.center[:,i] = (self.value[:,i]+self.value[:,i-1])/2.
        elif self.ftype == 'corner':
            for j in range(1,self.mesh.ny+1):
                for i in range(1,self.mesh.nx+1):
                    self.center[j,i] = (self.value[j,i]+self.value[j,i-1]+self.value[j-1,i]+self.value[j-1,i-1])/4.
        center.value = self.center 
        return center

    def compute_vedge(self):
        vedge = Fields_2d(self.mesh,ftype = 'vedge')
        self.vedge = np.zeros(self.mesh.center_x.shape)
        if self.ftype == 'center':
            for i in range(0,self.mesh.nx+1):
                self.vedge[:,i] = (self.value[:,i+1]+self.value[:,i])/2.
        elif self.ftype == 'hedge':
            for i in range(0,self.mesh.nx+1):
                for j in range(1,self.mesh.ny+2):
                    self.vedge[j,i] = (self.value[j-1,i]+self.value[j-1,i+1]+self.value[j,i]+self.value[j,i+1])/4.
        elif self.ftype == 'vedge':
            self.vedge = self.value
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        vedge.value = self.vedge
        return vedge

    def compute_hedge(self):
        hedge = Fields_2d(self.mesh,ftype = 'hedge')
        self.hedge = np.zeros(self.mesh.center_x.shape)
        if self.ftype == 'center':
            for j in range(0,self.mesh.ny+1):
                self.hedge[j,:] = (self.value[j,:]+self.value[j+1,:])/2.
        elif self.ftype == 'hedge':
            self.hedge = self.value
        elif self.ftype == 'vedge':
            for i in range(1,self.mesh.nx+2):
                for j in range(0,self.mesh.ny+1):
                    self.hedge[j,i] = (self.value[j,i]+self.value[j,i-1]+self.value[j+1,i]+self.value[j+1,i-1])/4.
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        hedge.value = self.hedge
        return hedge

    def compute_corner(self):
        corner = Fields_2d(self.mesh,ftype = 'corner')
        self.corner = np.zeros(self.mesh.center_x.shape)
        if self.ftype == 'center':
            for j in range(0,self.mesh.ny+1):
                for i in range(0,self.mesh.nx+1):
                    self.corner[j,i] = (self.value[j,i]+self.value[j,i+1]+self.value[j+1,i]+self.value[j+1,i+1])/4.
        elif self.ftype == 'hedge':
            for i in range(0,self.mesh.nx+1):
                self.corner[:,i] = (self.value[:,i+1]+self.value[:,i])/2.
        elif self.ftype == 'vedge':
            for j in range(0,self.mesh.ny+1):
                self.corner[j,:] = (self.value[j+1,:]+self.value[j,:])/2.
        elif self.ftype == 'corner':
            pass #currently there is no corner type
        corner.value = self.corner
        return corner

    def ddx_center(self):#d(self)/dx on CENTER!
        ddx_center = Fields_2d(self.mesh,ftype = 'center')
        if self.ftype == 'vedge':
            for i in range(1,self.mesh.nx+2):
                ddx_center.value[:,i] = (self.value[:,i]-self.value[:,i-1])/(self.mesh.dx)
        elif self.ftype == 'hedge':
            pass#todo
        elif self.ftype == 'center':
            pass#todo
        elif self.ftype == 'corner':
            pass #todo
        else:
            print "Fields not recognized when conducting ddx operation on ",self.name,". EXIT"
            sys.exit()
        return ddx_center

    def ddx_vedge(self):#d(self)/dx on VERTICAL EDGE
        ddx = Fields_2d(self.mesh,ftype = 'vedge')
        if self.ftype == 'vedge':
            pass #todo
        elif self.ftype == 'hedge':
            pass#todo
        elif self.ftype == 'center':
            for i in range(0,self.mesh.nx+1):
                ddx.value[:,i] = (self.value[:,i+1]-self.value[:,i])/(self.mesh.dx)
        elif self.ftype == 'corner':
            pass #todo
        else:
            print "Fields not recognized when conducting ddx operation on ",self.name,". EXIT"
            sys.exit()
        return ddx

    def ddy_center(self):#d(self)/dx
        ddy = Fields_2d(self.mesh,ftype = 'center')
        if self.ftype == 'vedge':
            pass#todo
        elif self.ftype == 'hedge':
            for j in range(1,self.mesh.ny+2):
                ddy.value[j,:] = (self.value[j,:]-self.value[j-1,:])/(self.mesh.dy)
        elif self.ftype == 'center':
            pass#todo
        elif self.ftype == 'corner':
            pass #todo
        else:
            print "Fields not recognized when conducting ddy operation on ",self.name,". EXIT"
            sys.exit()
        return ddy

    def ddy_hedge(self):#d(self)/dy on HORIZONTAL EDGE
        ddy = Fields_2d(self.mesh,ftype = 'hedge')
        if self.ftype == 'vedge':
            pass #todo
        elif self.ftype == 'hedge':
            pass#todo
        elif self.ftype == 'center':
            for j in range(0,self.mesh.ny+1):
                ddy.value[j,:] = (self.value[j+1,:]-self.value[j,:])/(self.mesh.dy)
        elif self.ftype == 'corner':
            pass #todo
        else:
            print "Fields not recognized when conducting ddx operation on ",self.name,". EXIT"
            sys.exit()
        return ddy

    def ddx2(self):#d2(self)/dx2 on self.ftype
        ddx2 = Fields_2d(self.mesh,ftype = self.ftype)
        for i in range(1,self.mesh.nx+1):
            ddx2.value[:,i] = (self.value[:,i+1]-2*self.value[:,i]+self.value[:,i-1])/(self.mesh.dx**2.)
        return ddx2
    def ddy2(self):#d2(self)/dy2 on self.ftype
        ddy2 = Fields_2d(self.mesh,ftype = self.ftype)
        for j in range(1,self.mesh.ny+1):
            ddy2.value[j,:] = (self.value[j+1,:]-2*self.value[j,:]+self.value[j-1,:])/(self.mesh.dy**2.)
        return ddy2

    def compute_residual(self,fields_2d):
        if self.ftype != fields_2d.ftype:
            print "cannot compute residual from two fields of different type!"
            sys.exit()
        residual  = 0
        for i in range(self.mesh.nx+2):
            for j in range(self.mesh.ny+2):
                residual = max(residual, abs(self.value[j,i]-fields_2d.value[j,i]))
        return residual


    
    def __add__(self,fields):
        #todo add a function to check if two fields can be added
        if (self.ftype != fields.ftype): 
            print "try to add up:",self.name,"and",fields.name
            print "cannot add two fields of different types up!"
            sys.exit()
        else:
        #def Fields_2d(self,mesh_2d, BC =None, IC='default',ftype='center', name='default_fields'):
            result = copy.deepcopy(self)
            result.value = self.value + fields.value
            return result

    def __sub__(self,fields):
        #todo add a function to check if two fields can be added
        if (self.ftype != fields.ftype): 
            print "try to add up:",self.name,"and",fields.name
            print "cannot add two fields of different types up!"
            sys.exit()
        else:
        #def Fields_2d(self,mesh_2d, BC =None, IC='default',ftype='center', name='default_fields'):
            result = copy.deepcopy(self)
            result.value = self.value - fields.value
            return result
    def __mul__(self,scalar):
        if type(scalar) != float:
            print "multiply operation must be taken between floats and fields"
            sys.exit()
        else:
            result = copy.deepcopy(self)
            #result = Fields_2d(self.mesh,self.BC,ftype=self.ftype)
            result.value = self.value * scalar
            return result
    def __rmul__(self,scalar):
        if type(scalar) != float:
            print "multiply operation must be taken between floats and fields"
            sys.exit()
        else:
            result = copy.deepcopy(self)
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
        #file.write("#"+self.name+'\n')
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

        


