#!/usr/bin/python
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from boundary_conditions import *
from initial_conditions import *
from poisson_solver import *
from multiphase import *
from fields_2d_velocity import *

meshx = Mesh_1d(0.,1.,0.025,0.025,'mesh_x')
meshy = Mesh_1d(0.,1.,0.025,0.025,'mesh_y')
mesh = Mesh_2d(meshx,meshy,'mesh_for_translation_test')
print "Mesh created.\n"
print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
mesh.write()
c=Fields_2d(mesh,ftype='center',IC=IC_translation_test,name='c')
[mx,my] = get_normal_young(c)#compute mx and my of each cell
alpha = get_alpha_field(c,mx,my)#compute alpha of each cell
draw_surface(c,mx,my,alpha,'translation_test_origin')#draw reconstructed circle

CFL = 0.25
u_inf = 1
max_iterations = 500000 
u = U(mesh,IC=IC_translation_test_u_1)
v = V(mesh,IC=IC_translation_test_v_1)
t=0.0
u_max = max(np.amax(np.fabs(u.value)),np.amax(np.fabs(v.value)))
#dt=CFL*mesh.min_delta/u_max
dt = 0.005
CFL = u_max*dt/mesh.min_delta
frame = 0
output_interval= 200 #Interval that fields are written
x_first = True #sweep in x direction first


c.write(str(t))
print "dt=",dt
print "CFL=",CFL


#main loop
for n in range(max_iterations):
    frame = frame +1
    t= t+dt
    if t >= 1.0 and t<2.0:
        u.value = fixed_value(mesh,0)
        v.value = fixed_value(mesh,-0.2)
    elif t >= 2.0 and t<3.0:
        u.value = fixed_value(mesh,-0.8)
        v.value = fixed_value(mesh,0.4)
    elif t >= 3.0 and t<4.0:
        u.value = fixed_value(mesh,0.4)
        v.value = fixed_value(mesh,-0.4)
    elif t >= 4.0 and t<5.0:
        u.value = fixed_value(mesh,-0.2)
        v.value = fixed_value(mesh,0.0)
    elif t >= 5.0 and t<6.0:
        u.value = fixed_value(mesh,0.2)
        v.value = fixed_value(mesh,0.8)
    if t > 6.0:
        print "loop ending time:",t
        break
    print "iteration:",n+1
    if x_first == True:
        x_first = False
        #sweep in x direction
        fx_c = u.uphi_x_vedge_c(c,mx,my,alpha,dt)
        fx_c.name = 'fx_c'
        fx_c_x = fx_c.ddx_center()
        c.value = (-1)*fx_c_x.value*mesh.dx+c.value

        #reconstruct
        [mx,my] = get_normal_young(c)#compute mx and my of each cell
        alpha = get_alpha_field(c,mx,my)#compute alpha of each cell

        #sweep in y direction
        fy_c = v.vphi_y_hedge_c(c,mx,my,alpha,dt)
        fy_c.name = 'fy_c'
        fy_c_y = fy_c.ddy_center()
        c.value = (-1)*fy_c_y.value*mesh.dy+c.value
        
        #reconstruct
        [mx,my] = get_normal_young(c)#compute mx and my of each cell
        alpha = get_alpha_field(c,mx,my)#compute alpha of each cell
    else:
        x_first = True
        #sweep in y direction
        fy_c = v.vphi_y_hedge_c(c,mx,my,alpha,dt)
        fy_c.name = 'fy_c'
        fy_c_y = fy_c.ddy_center()
        c.value = (-1)*fy_c_y.value*mesh.dy+c.value
        
        #reconstruct
        [mx,my] = get_normal_young(c)#compute mx and my of each cell
        alpha = get_alpha_field(c,mx,my)#compute alpha of each cell

        #sweep in x direction
        fx_c = u.uphi_x_vedge_c(c,mx,my,alpha,dt)
        fx_c.name = 'fx_c'
        fx_c_x = fx_c.ddx_center()
        c.value = (-1)*fx_c_x.value*mesh.dx+c.value
        #reconstruct
        [mx,my] = get_normal_young(c)#compute mx and my of each cell
        alpha = get_alpha_field(c,mx,my)#compute alpha of each cell

    if frame == output_interval:
        fx_c.write(str(t))
        fy_c.write(str(t))
        c.write(str(t))
        draw_surface(c,mx,my,alpha,'circle_at_t='+str(t))#draw reconstructed circle
        frame = 0


