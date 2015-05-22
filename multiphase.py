#!/usr/bin/python
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from boundary_conditions import *
from fields_2d_velocity import U,V
from initial_conditions import *
from poisson_solver import *
import numpy as np
import math
import matplotlib.pyplot as plt
import pdb

#meshx = Mesh_1d(1.,4.,1,1,'mesh_x')
#meshy = Mesh_1d(1.,3.,1,1,'mesh_y')
#mesh = Mesh_2d(meshx,meshy,'mesh')
#print "Mesh created.\n"
#print "Mesh size: (",mesh.nx,', ',mesh.ny,")"
#mesh.write()
#c=Fields_2d(mesh,ftype='center',name='c')
#c.value = np.array([[1,1,1,1,1],[1./3.,11./12.,1,1,1],[0,1./12,2./3,1,1],[0,0,0,1./3,11./12]])
##above are for testing young's method



def get_normal_young(c):#return mx my from c field
    c_x = c.ddx_vedge()
    mx_corner = -1. * c_x.compute_corner()
    mx = mx_corner.compute_center()
    c_y=c.ddy_hedge()
    my_corner=-1.*c_y.compute_corner()
    my=my_corner.compute_center()
    #alpha = Fields_2d(mesh,ftype='center',name = 'alpha')
    #k=-1*mx.value/my.value
    return mx,my

def get_alpha_field(c,mx,my):#c,mx,my are Fields_2d objects
    alpha = Fields_2d(c.mesh,ftype = 'center',name='alpha')
    for j in range(1,c.mesh.ny+1):
        for i in range(1,c.mesh.nx+1):
            if (1 - c.value[j,i]) >0.01 and (c.value[j,i]-0) > 0.01:#surface cell
                if mx.value[j,i]>0  and my.value[j,i] >0:
                    alpha.value[j,i] = compute_alpha(c.mesh,mx.value[j,i],my.value[j,i],c.value[j,i])
                elif mx.value[j,i]>0 and my.value[j,i] <0:
                    alpha_bar = compute_alpha(c.mesh,mx.value[j,i],-1*my.value[j,i],c.value[j,i])
                    alpha.value[j,i] = alpha_bar - 2*(-1)*my.value[j,i]*c.mesh.dy/2.
                elif mx.value[j,i]<0 and my.value[j,i]>0:
                    alpha_bar = compute_alpha(c.mesh,-1*mx.value[j,i],my.value[j,i],c.value[j,i])
                    alpha.value[j,i] = alpha_bar - 2*(-1)*mx.value[j,i]*c.mesh.dx/2.
                elif mx.value[j,i]<0 and my.value[j,i]<0:
                    alpha_bar = compute_alpha(c.mesh,-1*mx.value[j,i],-1*my.value[j,i],c.value[j,i])
                    alpha.value[j,i] = alpha_bar - 2*(-1)*mx.value[j,i]*c.mesh.dx/2. - 2*(-1)*my.value[j,i]*c.mesh.dy/2.
                #else:
                #    print "There may be something wrong with input mx and my"
                #    sys.exit()
    return alpha


#[mx,my] = reconstruct_surface_young(c)

def area(mesh,alpha,mx,my):
    area = 1./(2*mx*my)*(alpha**2-max(0,alpha-mx*mesh.dx)**2-max(0,alpha-my*mesh.dy)**2)
    if area > mesh.dx*mesh.dy*1.1 or area<0:
        print "area is out of range: (0,dx*dy). Something must be wrong"
        sys.exit()
    return area

#finish testing
def compute_alpha(mesh,mxi,myi,ci):#mxi and myi must be positive
    if mxi<0 or myi <0 :
        print "imput mxi,myi and ci must be positive"
        sys.exit()
    #alpha = (mxi*mesh.dx+myi*mesh.dy)/2. #initial guess
    max_iteration = 500
    convergence = 0.0001
    #A = area(mesh,alpha,mxi,myi)
    ci = ci*mesh.dx*mesh.dy
    print "target area: ",ci
    print "mx,my:",mxi,myi
    #for n in range(max_iteration):#find two alpha,one larger than root, one less than root
        #if A > mesh.dx*mesh.dy*1.1:
        #    pdb.set_trace()
        #print alpha
        #print "area:",A
        ##pdb.set_trace()
        #if A>ci:
        #    alpha2 = alpha
        #    alpha = 0.5*alpha
        #elif A<ci:
        #    alpha2 = alpha
        #    alpha = 1.5*alpha
        #A2 = A#(alpha2 -> A2)
      
        #A = area(mesh,alpha,mxi,myi)
        #if (A2-ci)*(A-ci)<0:
        #    break
    alpha = 0.
    alpha2 = (mxi*mesh.dx+myi*mesh.dy)
    A = area(mesh,alpha,mxi,myi)

    #pdb.set_trace()
    #print A2,A,alpha2,alpha
    #Now one of A and A2 is above ci and one is below
    alpha3 = (alpha2+alpha)/2.
    for n in range(max_iteration):
        #if abs(ci-0.0594135802469) < 0.001:
        #    pdb.set_trace()
        print "alpha3:",alpha3
        print "alpha :",alpha
        print "alpha2:",alpha2
        A3 = area(mesh,alpha3,mxi,myi)
        if (A3-ci)*(A-ci)<0:
            alpha2 = alpha3
            alpha3 = (alpha3+alpha)/2.
        else:
            alpha = alpha3
            alpha3 = (alpha3+alpha2)/2.
            A = area(mesh,alpha,mxi,myi)
        residual = abs(A3-ci)/ci
        print "iteration:",n+1
        print "residual:",residual
        #if math.isnan(residual):
        #    pdb.set_trace()
        
        if residual<convergence:
            return alpha
    #pdb.set_trace()


    print "cannot get a proper alpha, in function compute_alpha()"
    sys.exit()

def draw_surface(c,mx,my,alpha,name=''):
    plt.figure() 
    for j in range(1,c.mesh.ny+1):
        for i in range(1,c.mesh.nx+1):
            if (1 - c.value[j,i]) >0.01 and (c.value[j,i]-0) > 0.01:#surface cell
                popout =[]
                #x,y in local coordinates
                x = np.linspace(0,c.mesh.dx,20)
                y = -1*mx.value[j,i]/my.value[j,i]*x+alpha.value[j,i]/my.value[j,i]
                for k in range(x.size):
                    if y[k]>c.mesh.dy or y[k]<0:
                        popout.append(k)
                x = np.delete(x,popout)
                #transfer to global coordinates
                x = x+c.mesh.vedge_x[j,i-1]
                y = -1*mx.value[j,i]/my.value[j,i]*x+(alpha.value[j,i]+mx.value[j,i]*c.mesh.vedge_x[j,i-1]+my.value[j,i]*c.mesh.hedge_y[j-1,i])/my.value[j,i]
                plt.plot(x,y,'-')
    #plot mesh grid
    for j in range(0,c.mesh.ny+1):
        x = np.linspace(c.mesh.x0,c.mesh.xn,50)
        y = np.zeros(x.size)+c.mesh.hedge_y[j,0]
        plt.plot(x,y,'b-')
    for i in range(0,c.mesh.nx+1):
        y = np.linspace(c.mesh.y0,c.mesh.yn,50)
        x = np.zeros(y.size)+c.mesh.vedge_x[0,i]
        plt.plot(x,y,'b-')

    ax = plt.gca()
    #ax.grid(b=True, which='major', color='b', linestyle='-')
    #plt.minorticks_on()
    #ax.grid(b=True, which='minor', color='r', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.xlim((2,4))
    #plt.ylim((6,9))
    ax.set_aspect('equal')
    if not ('myplot' in os.listdir('./')):
        os.mkdir('./myplot')
    plt.savefig('./myplot/surface_'+name+'.png', bbox_inches='tight')
    plt.close()






