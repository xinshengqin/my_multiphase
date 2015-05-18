#!/usr/bin/python
from mesh import mesh_1d,mesh_2d
from fields import fields_2d 
from boundary_conditions import *


mesh1 = mesh_1d(0.,6.,2.,2.)
#mesh1.plot_edgeVsIndex()
#mesh1.plot_centerVsIndex()

mesh2 = mesh_1d(0,0.1,0.05,0.05)
#mesh2.plot_edgeVsIndex()
#mesh2.plot_centerVsIndex()


mesh3 = mesh_2d(mesh1,mesh2,'mesh3')
mesh3.write()

u = fields_2d(mesh3,BC_fixedValue_1)
u.write('before_BC')
u.apply_BC()
u.write('1.1')


