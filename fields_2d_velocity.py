#!/usr/bin/python
import numpy as np
import sys
import os
from fields import Fields_2d
from multiphase import compute_flux_c_lr,compute_flux_c_rl,compute_flux_c_lu,compute_flux_c_ul

class U(Fields_2d):
    def __init__(self,mesh_2d, BC = None, IC = 'default', name = 'U'):
        Fields_2d.__init__(self,mesh_2d,BC,IC, ftype = 'vedge',name = 'U')

    def uphi_x_vedge(self,phi):#\frac{\partial u*phi}{\partial x} on vertical edge
        #phi is prefered to be defined on vertical edge
        uphix = Fields_2d(self.mesh,ftype = 'vedge')
        self.compute_center()
        phi.compute_center()
        for i in range(0,self.mesh.nx+1):
            uphix.value[:,i] = 1./self.mesh.dx * (self.center[:,i+1]*phi.center[:,i+1]-self.center[:,i]*phi.center[:,i])
        return uphix

    def uphi_x_hedge(self,phi):#\frac{\partial u*phi}{\partial x} on horizontal edge
        #phi prefered to be defined on horizontal edge
        uphix = Fields_2d(self.mesh,ftype = 'hedge')
        self.compute_corner()
        phi.compute_corner()
        for i in range(1,self.mesh.nx+2):
            uphix.value[:,i] = 1./self.mesh.dx * (self.corner[:,i]*phi.corner[:,i]-self.corner[:,i-1]*phi.corner[:,i-1])
        return uphix

    def uphi_x_vedge_c(self,c,mx,my,alpha,dt):
        fx_c = Fields_2d(self.mesh,ftype = 'vedge')#flux of c in x direction
        for j in range(1,self.mesh.ny+1):
            for i in range(1,self.mesh.nx+1):
                if self.value[j,i]>0:
                    if (c.value[j,i]-0.)>0.001 and (1-c.value[j,i])>0.001:#surface cell
                    #compute flux of c from left to right
                    #use cell at left side of this edge to compute flux
                        fx_c.value[j,i] = compute_flux_c_lr(mx.value[j,i],my.value[j,i],alpha.value[j,i],self.value[j,i]*dt,self.mesh.dx,self.mesh.dy)
                    elif (c.value[j,i]-0.)<0.001:#empty cell 
                        fx_c.value[j,i] = 0
                    elif (1-c.value[j,i])<0.001:#full cell
                        fx_c.value[j,i] = (self.value[j,i]*dt*self.mesh.dy)/(self.mesh.dx*self.mesh.dy)

                elif self.value[j,i]<0:#u[i+1/2,j]<0
                #compute flux of c from right to left
                #use cell at right side of this edge to compute flux
                    if (c.value[j,i+1]-0.)>0.001 and (1-c.value[j,i+1])>0.001:#surface cell
                        fx_c.value[j,i] = compute_flux_c_rl(mx.value[j,i+1],my.value[j,i+1],alpha.value[j,i+1],self.value[j,i]*dt,self.mesh.dx,self.mesh.dy)
                    elif (c.value[j,i+1]-0.)<0.001:#empty cell 
                        fx_c.value[j,i] = 0
                    elif (1-c.value[j,i+1])<0.001:#full cell
                        fx_c.value[j,i] = (self.value[j,i]*dt*self.mesh.dy)/(self.mesh.dx*self.mesh.dy)
        return fx_c







class V(Fields_2d):
    def __init__(self,mesh_2d, BC = None, IC = 'default', name = 'V'):
        Fields_2d.__init__(self,mesh_2d,BC, IC, ftype = 'hedge',name = 'V')
        self.BC = BC

    def vphi_y_vedge(self,phi):#\frac{\partial v*phi}{\partial y} on vertical edge
        vphiy = Fields_2d(self.mesh,ftype='vedge')
        self.compute_corner()
        phi.compute_corner()
        for j in range(1,self.mesh.ny+2):
            vphiy.value[j,:] = 1./self.mesh.dy * (self.corner[j,:]*phi.corner[j,:]-self.corner[j-1,:]*phi.corner[j-1,:])
        return vphiy


    def vphi_y_hedge(self,phi):#\frac{\partial v*phi}{\partial y} on vertical edge
        vphiy = Fields_2d(self.mesh,ftype='hedge')
        self.compute_center()
        phi.compute_center()
        for j in range(0,self.mesh.ny+1):
            vphiy.value[j,:] = 1./self.mesh.dy * (self.center[j+1,:]*phi.center[j+1,:]-self.center[j,:]*phi.center[j,:])
        return vphiy

    def vphi_y_hedge_c(self,c,mx,my,alpha,dt):
        fy_c = Fields_2d(self.mesh,ftype = 'hedge')#flux of c in y direction
        for j in range(1,self.mesh.ny+1):
            for i in range(1,self.mesh.nx+1):
                if self.value[j,i]>0:
                    if (c.value[j,i]-0.)>0.001 and (1-c.value[j,i])>0.001:#surface cell
                    #compute flux of c from lower to upper 
                    #use cell at lower side of this edge to compute flux
                        fy_c.value[j,i] = compute_flux_c_lu(mx.value[j,i],my.value[j,i],alpha.value[j,i],self.value[j,i]*dt,self.mesh.dx,self.mesh.dy)
                    elif (c.value[j,i]-0.)<0.001:#empty cell 
                        fy_c.value[j,i] = 0
                    elif (1-c.value[j,i])<0.001:#full cell
                        fy_c.value[j,i] = (self.value[j,i]*dt*self.mesh.dx)/(self.mesh.dx*self.mesh.dy)

                elif self.value[j,i]<0:#v[i,j+1/2]<0
                #compute flux of c from upper to lower
                #use cell at upper side of this edge to compute flux
                    if (c.value[j+1,i]-0.)>0.001 and (1-c.value[j+1,i])>0.001:#surface cell
                        fy_c.value[j,i] = compute_flux_c_ul(mx.value[j+1,i],my.value[j+1,i],alpha.value[j+1,i],self.value[j,i]*dt,self.mesh.dx,self.mesh.dy)
                    elif (c.value[j+1,i]-0.)<0.001:#empty cell 
                        fy_c.value[j,i] = 0
                    elif (1-c.value[j+1,i])<0.001:#full cell
                        fy_c.value[j,i] = (self.value[j,i]*dt*self.mesh.dx)/(self.mesh.dx*self.mesh.dy)
        return fy_c
