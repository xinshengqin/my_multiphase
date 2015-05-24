#!/usr/bin/python
from mesh import Mesh_1d,Mesh_2d
from fields import Fields_2d 
from boundary_conditions import *
from initial_conditions import *
from poisson_solver import *
import numpy as np
import math
import matplotlib.pyplot as plt
import pdb

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
            if (1 - c.value[j,i]) >0.001 and (c.value[j,i]-0) > 0.001:#surface cell
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

def area(alpha,mx,my,dx,dy):#alpha,mx,my are floats 
    #mx,my must be positive
    if mx < 0 or my < 0:
        print "In function area(), input mx, my must be positive"
        sys.exit()
    area = 1./(2*mx*my)*(alpha**2-max(0,alpha-mx*dx)**2-max(0,alpha-my*dy)**2)
    #if area > dx*dy*1.0:
    if alpha > mx*dx+my*dy: 
        return dx*dy
    if area < 0.0:
        print "area is negative. Something must be wrong"
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
    #print "target area: ",ci
    #print "mx,my:",mxi,myi
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
    #A = area(mesh,alpha,mxi,myi)
    A = area(alpha,mxi,myi,mesh.dx,mesh.dy)

    #pdb.set_trace()
    #print A2,A,alpha2,alpha
    #Now one of A and A2 is above ci and one is below
    alpha3 = (alpha2+alpha)/2.
    for n in range(max_iteration):
        #if abs(ci-0.0594135802469) < 0.001:
        #    pdb.set_trace()
        #print "alpha3:",alpha3
        #print "alpha :",alpha
        #print "alpha2:",alpha2
        A3 = area(alpha3,mxi,myi,mesh.dx,mesh.dy)
        if (A3-ci)*(A-ci)<0:
            alpha2 = alpha3
            alpha3 = (alpha3+alpha)/2.
        else:
            alpha = alpha3
            alpha3 = (alpha3+alpha2)/2.
            A = area(alpha,mxi,myi,mesh.dx,mesh.dy)
        residual = abs(A3-ci)/ci
        #print "iteration:",n+1
        #print "residual:",residual
        #if math.isnan(residual):
        #    pdb.set_trace()
        
        if residual<convergence:
            return alpha
    #pdb.set_trace()


    print "cannot get a proper alpha, in function compute_alpha()"
    sys.exit()

#finish testing
def compute_flux_c_lr(mxi,myi,alphai,uidt,dx,dy):#uidt is ui*dt
    #pdb.set_trace()
    flux = 0.
    if uidt<0:
        print "ui must be positive in compute_flux_c_lr()"
        sys.exit()

    if mxi>0 and myi>0:
        flux = max(0,(area(alphai,mxi,myi,dx,dy)-area(alphai,mxi,myi,dx-uidt,dy))/(dx*dy))
    elif mxi>0 and myi<0:
        alphai = alphai-2*myi*dy/2.
        myi = -1*myi
        flux = max(0,(area(alphai,mxi,myi,dx,dy)-area(alphai,mxi,myi,dx-uidt,dy))/(dx*dy))
    elif mxi<0 and myi>0:
        alphai = alphai - 2*mxi*dx/2.
        mxi = -1*mxi
        flux = max(0,area(alphai,mxi,myi,uidt,dy)/(dx*dy))
    elif mxi<0 and myi<0:
        alphai = alphai - 2*mxi*dx/2. - 2*myi*dy/2.
        mxi = -1*mxi
        myi = -1*myi
        #print flux
        flux = max(0,area(alphai,mxi,myi,uidt,dy)/(dx*dy))

    if flux>1.:
        print "flux cannot be larger than 1 from compute_flux_c_lr()"
        sys.exit()
    return flux

#finish testing
def compute_flux_c_rl(mxi,myi,alphai,uidt,dx,dy):
    flux = 0.
    if uidt>0:
        print "ui must be negative in compute_flux_c_rl()"
        sys.exit()
    if mxi>0 and myi>0:
        flux = max(0,area(alphai,mxi,myi,-1*uidt,dy)/(dx*dy))
    elif mxi>0 and myi<0:
        alphai = alphai - 2*myi*dy/2.
        myi = -1*myi
        flux = max(0,area(alphai,mxi,myi,-1*uidt,dy)/(dx*dy))
    elif mxi<0 and myi>0:
        alphai = alphai - 2*mxi*dx/2.
        mxi = -1*mxi
        flux = max(0,(area(alphai,mxi,myi,dx,dy)-area(alphai,mxi,myi,dx-(-1*uidt),dy))/(dx*dy))
        #a1 = area(alphai,mxi,myi,dx,dy)
        #a2 = area(alphai,mxi,myi,dx-(-1*uidt),dy)
        #pdb.set_trace()
    elif mxi<0 and myi<0:
        alphai = alphai - 2*mxi*dx/2. - 2*myi*dy/2.
        mxi = -1*mxi
        myi = -1*myi
        flux = max(0,(area(alphai,mxi,myi,dx,dy)-area(alphai,mxi,myi,dx-(-1*uidt),dy))/(dx*dy))
    if flux>1:
        print "flux cannot be less than -1 from compute_flux_c_rl()"
        sys.exit()
    return flux*(-1)

#compute flux of c from lower to upper 
def compute_flux_c_lu(mxi,myi,alphai,vidt,dx,dy):
    flux = 0.
    if vidt<0:
        print "vi must be positive in compute_flux_c_lu()"
        sys.exit()
    #rotate the system
    alpha_n = alphai - mxi*dx
    mx_n = myi
    my_n = -1*mxi
    dy_n = dx
    dx_n = dy
    flux = compute_flux_c_lr(mx_n,my_n,alpha_n,vidt,dx_n,dy_n)
    return flux


def compute_flux_c_ul(mxi,myi,alphai,uidt,dx,dy):
    flux = 0.
    if uidt>0:
        print "vi must be negative in compute_flux_c_ul()"
        sys.exit()
    #rotate the system
    alpha_n = alphai - mxi*dx
    mx_n = myi
    my_n = -1*mxi
    dy_n = dx
    dx_n = dy
    flux = compute_flux_c_rl(mx_n,my_n,alpha_n,uidt,dx_n,dy_n)
    return flux

def draw_surface(c,mx,my,alpha,name=''):
    plt.figure() 
    resolution = 20
    for j in range(1,c.mesh.ny+1):
        for i in range(1,c.mesh.nx+1):
            if (1 - c.value[j,i]) >0.001 and (c.value[j,i]-0) > 0.001:#surface cell
                popout =[]
                #x,y in local coordinates
                x = np.linspace(0,c.mesh.dx,resolution)
                y = -1*mx.value[j,i]/my.value[j,i]*x+alpha.value[j,i]/my.value[j,i]
                for k in range(x.size):
                    if y[k]>c.mesh.dy or y[k]<0:
                        popout.append(k)
                x = np.delete(x,popout)
                if x.size<4:#line is almost vertical
                    popout = []
                    y = np.linspace(0,c.mesh.dy,resolution)
                    x = -1*my.value[j,i]/mx.value[j,i]*y+alpha.value[j,i]/mx.value[j,i]
                    for k in range(y.size):
                        if x[k]>c.mesh.dx or x[k]<0:
                            popout.append(k)
                    x = np.delete(x,popout)#delete points out of this cell
                #transfer to global coordinates
                x = x+c.mesh.vedge_x[j,i-1]
                y = -1*mx.value[j,i]/my.value[j,i]*x+(alpha.value[j,i]+mx.value[j,i]*c.mesh.vedge_x[j,i-1]+my.value[j,i]*c.mesh.hedge_y[j-1,i])/my.value[j,i]
                plt.plot(x,y,'-',linewidth=0.5)
                    
    #plot mesh grid
    for j in range(0,c.mesh.ny+1):
        x = np.linspace(c.mesh.x0,c.mesh.xn,50)
        y = np.zeros(x.size)+c.mesh.hedge_y[j,0]
        plt.plot(x,y,'b-',linewidth=0.1)
    for i in range(0,c.mesh.nx+1):
        y = np.linspace(c.mesh.y0,c.mesh.yn,50)
        x = np.zeros(y.size)+c.mesh.vedge_x[0,i]
        plt.plot(x,y,'b-',linewidth=0.1)

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
    plt.savefig('./myplot/surface_'+name+'.png', bbox_inches='tight',dpi=200)
    plt.close()






