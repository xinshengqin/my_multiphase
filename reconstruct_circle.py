# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from boundary_conditions import *
from initial_conditions import *
from poisson_solver import *
from multiphase import *

meshx = Mesh_1d(0.,10.,0.25,0.25,'mesh_x')
meshy = Mesh_1d(0.,10.,0.25,0.25,'mesh_y')
mesh = Mesh_2d(meshx,meshy,'mesh_for_reconstruct_circle')
print "Mesh created.\n"
print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
mesh.write()
c=Fields_2d(mesh,ftype='center',IC=IC_C_circle,name='c')
c.write()
[mx,my] = get_normal_young(c)#compute mx and my of each cell
alpha = get_alpha_field(c,mx,my)#compute alpha of each cell
draw_surface(c,mx,my,alpha,'circle')#draw reconstructed circl
