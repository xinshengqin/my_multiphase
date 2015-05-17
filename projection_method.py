#!/usr/bin/python
from mesh import mesh_1d,mesh_2d
from fields import fields_2d 
from boundary_conditions import *


#create 2d mesh
meshx = mesh_1d(0.,10.,0.1,0.1,'mesh_x')
meshy = mesh_1d(0.,4.,0.1,0.1,'mesh_y')
meshx.plot_edgeVsIndex()
meshx.plot_centerVsIndex()
meshy.plot_edgeVsIndex()
meshy.plot_centerVsIndex()
mesh = mesh_2d(meshx,meshy,'mesh')
print "Mesh created.\n"
print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
mesh.write()

CFL = 0.8
u_inf = 1
max_iterations = 10000
#todo defined a IC function for this problem
#todo modify BC
u = fields_2d(mesh,BC_fixedValue_1,IC)
v = fields_2d(mesh,BC_fixedValue_1,IC)
p = fields_2d(mesh,BC_fixedValue_1,IC)
fields_list = [u,v,p]

n_steps=10000 #Interval that u,v and p are write
Re=100.0   #Reynolds number
nu=1e-6
rho=1000.0 
t=0.0
dt=CFL*mesh.min_delta/u_inf

#write all fields into time directory
for i,item in enumerate(fields_list):
    item.write(str(t))

flux_ux = fields_2d(mesh,BC_empty,IC='default')
flux_uy = fields_2d(mesh,BC_empty,IC='default')
flux_vy = fields_2d(mesh,BC_empty,IC='default')
flux_vx = fields_2d(mesh,BC_empty,IC='default')

#main loop
for n in range(max_iterations):
    print "iteration: ",n,"\n"
    u_old = u
    v_old = v
    for i,item in enumerate(fields_list):
        item.applyBC()


    # Step 1: Update velocity to intermediate step
    # Calculate fluxes at each boundary
    for i in range(2,mesh.nx+4):
        for j in range(0,mesh.ny+4): #j = 0,1,...mesh.ny+3
            flux_ux.value[i,j] = (u.value[i-1,j]+u.value[i,j])**2 / 4.

    for i in range(0,mesh.nx+2):
        for j in range(0,mesh.ny+2):
            flux_uy.value[i,j] = (u.value[i,j]+u.value[i,j+1]) * (v.value[i+1,j]+v.value[i,j]) / 4.

    for i in range(0,mesh.nx+4):
        for j in range(0,mesh.ny+2):
            flux_vy.value[i,j] = (v.value[i,j+1] + v.value[i,j])**2 / 4.

    for i in range(0,mesh.nx+2):
        for j in range(0,mesh.ny+2): #j = 0,1,...mesh.ny+3
            flux_vx.value[i,j] = (u.value[i,j+1]+u.value[i,j])*(v.value[i,j]+v.value[i+1,j]) / 4. 

    for i in range(2,mesh.nx+2):
        for j in range(2,mesh.ny+2):
            # Advective terms
            #in RU(i,j)
            uu_x = (flux_ux.value[i+1,j]-flux_ux.value[i,j])/(mesh.center_x[i+1,j]-mesh.center_x[i,j]) #\frac{\partial uu}{\partial x}
            uv_y = (flux_uy.value[i,j]-flux_uy.value[i,j-1])/(mesh.vedge_x[i,j]-mesh.vedge_x[i,j-1])
            
            #in RV(i,j)
            #uv_x = (Flux_vx(i,j)-Flux_vx(i-1,j))/dx
            uv_x = (flux_vx.value[i,j]-flux_vx.value[i-1,j])/(mesh.vedge_x[i,j]-mesh.vedge[i-1,j])
            vv_y = (flux_vy.value[i,j+1]-flux_vy.value[i,j])/(mesh.center_y[i,j+1]-mesh.center_y[i,j])
            
            # Diffusive terms
            u_x_ip1 = (u.value[i+1,j]-u.value[i,j])/(mesh.vedge_x[i+1]-mesh.vedge_x[i,j])
            u_x_i = (u.value[i,j]-u.value[i-1,j])/(mesh.vedge_x[i]-mesh.vedge_x[i-1,j])
            u_xx = 1./(mesh.center_x[i+1,j]-mesh.center_x[i,j])*(u_x_ip1 - u_x_i)


            u_yy = F_center(j)*(F_edge(j)*(u(i,j+1)-u(i,j))-F_edge(j-1)*(u(i,j)-u(i,j-1)))/(dzeta**2)
            v_xx = (v(i+1,j)-2*v(i,j)+v(i-1,j))/(dx**2)
            v_yy = F_edge(j)*(F_center(j+1)*(v(i,j+1)-v(i,j))-F_center(j)*(v(i,j)-v(i,j-1)))/(dzeta**2)
            
            ! Update to u* and v* values
            u_star(i,j) = u(i,j)+dt*(-(uu_x+uv_y)+1/Re*(u_xx+u_yy))
            v_star(i,j) = v(i,j)+dt*(-(uv_x+vv_y)+1/Re*(v_xx+v_yy))
