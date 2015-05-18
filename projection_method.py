#!/usr/bin/python
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from fields_2d_velocity import U,V
from boundary_conditions import *


#create 2d mesh
meshx = Mesh_1d(0.,10.,0.1,0.1,'mesh_x')
meshy = Mesh_1d(0.,4.,0.1,0.1,'mesh_y')
#meshx.plot_edgeVsIndex()
#meshx.plot_centerVsIndex()
#meshy.plot_edgeVsIndex()
#meshy.plot_centerVsIndex()
mesh = Mesh_2d(meshx,meshy,'mesh')
print "Mesh created.\n"
print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
mesh.write()

CFL = 0.8
u_inf = 1
max_iterations = 10000
#todo defined a IC function for this problem
#todo modify BC
u = U(mesh,BC_fixedValue_1,IC)
v = V(mesh,BC_fixedValue_1,IC)
p = Fields_2d(mesh,BC_fixedValue_1,IC)
fields_list = [u,v,p]

n_steps=10000 #Interval that u,v and p are write
Re=100.0   #Reynolds number
#nu=1e-6
mu = 1e-3
rho=1000.0 
t=0.0
dt=CFL*mesh.min_delta/u_inf

#write all fields into time directory
for i,item in enumerate(fields_list):
    item.write(str(t))

u_star = U(mesh)
v_star = V(mesh)

#main loop
for n in range(max_iterations):
    print "iteration: ",n,"\n"
    u_old = u
    v_old = v
    for i,item in enumerate(fields_list):
        item.applyBC()


    # Step 1: Update velocity to intermediate step
    # Calculate fluxes at each boundary
    Ax = u.uphi_x_vedge(u)+v.vphi_y_vedge(u)
    Dx = mu*(u.ddx2()+u.ddy2())
    Ay = u.uphi_x_hedge(v)+v.vphi_y_hedge(v)
    Dy = mu*(v.ddx2()+v.ddy2())
    #todo
    fbx = Fields_2d(mesh,ftype='vedge')#currently its zero
    fby = Fields_2d(mesh,ftype='hedge')#currently its zero
            
    u_star.value = u.value + dt*(-1.*Ax+fbx+1./rho*(Dx))
    v_star.value = v.value + dt*(-1.*Ay+fby+1./rho*(Dy))

    break
    
            



