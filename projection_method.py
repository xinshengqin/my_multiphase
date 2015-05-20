#!/usr/bin/python
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from boundary_conditions import *
from fields_2d_velocity import U,V
from initial_conditions import *
from poisson_solver import *
import numpy as np

para={}
para['rho']= 1.0
para['max_teration_poission'] = 10000
para['tolerance_poisson'] = 1e-4
para['w'] = 1.6


#create 2d mesh
meshx = Mesh_1d(0.,10.,0.1,0.1,'mesh_x')
meshy = Mesh_1d(0.,4,0.1,0.1,'mesh_y')
meshx.plot_edgeVsIndex()
meshx.plot_centerVsIndex()
meshy.plot_edgeVsIndex()
meshy.plot_centerVsIndex()
mesh = Mesh_2d(meshx,meshy,'mesh')
print "Mesh created.\n"
print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
mesh.write()

CFL = 0.8
u_inf = 1
max_iterations = 10000
u = U(mesh,BC_flatplate_u,IC_flatplate_u)
v = V(mesh,BC_flatplate_u,IC_flatplate_v)
p = Fields_2d(mesh,BC_p,IC_flatplate_p,name='p')
fields_list = [u,v,p]

output_interval=1 #Interval that u,v and p are write
tolerance = 0.001
Re=100.0   #Reynolds number
#nu=1e-6
mu = 0.05
rho=para['rho']
t=0.0
dt=CFL*meshx.min_delta/u_inf/Re
para['dt'] = dt

#write all fields into time directory
for i,item in enumerate(fields_list):
    item.write(str(t))

u_star = U(mesh)
v_star = V(mesh)
frame = 0
residual = 0

#main loop
for n in range(max_iterations):
    print "iteration: ",n,"\n"
    t= t+dt
    frame = frame + 1
    print "time = ",t,"\n" 
    u_old = u
    v_old = v
    for i,item in enumerate(fields_list):
        item.applyBC()

    # Step 1: Update velocity to intermediate step
    Ax = u.uphi_x_vedge(u)+v.vphi_y_vedge(u)
    Dx = mu*(u.ddx2()+u.ddy2())
    Ay = u.uphi_x_hedge(v)+v.vphi_y_hedge(v)
    Dy = mu*(v.ddx2()+v.ddy2())
    fbx = Fields_2d(mesh,ftype='vedge')#currently its zero
    fby = Fields_2d(mesh,ftype='hedge')#currently its zero
            
    u_star = u + (dt*(-1.*Ax+fbx+1./rho*(Dx)))
    v_star = v + (dt*(-1.*Ay+fby+1./rho*(Dy)))

    # Step 2: Solve projection poisson problem
    u_star.applyBC()
    v_star.applyBC()
    [status,p]=solve_poisson(mesh,u_star,v_star,p,para)


    # Step 3: Update velocities to end time
    u = u_star - (dt/rho/mesh.dx)*p.ddx_vedge()
    v = v_star - (dt/rho/mesh.dx)*p.ddy_hedge()
    # Step 4: Check convergence
    residual =  max(residual,u.compute_residual(u_old),v.compute_residual(v_old))
    print "residual = ",residual
    if residual < tolerance:
        print "Convergence reached, R = ",residual
        for i,item in enumerate(fields_list):
            item.write(str(t))
        break
    if frame == output_interval:
        for i,item in enumerate(fields_list):
            item.write(str(t))
        frame = 0
    
            



