#!/usr/bin/python
import sys

def solve_poisson(mesh,u,v,p,para):
    Q = (para['rho']/para['dt']) * (u.ddx_center()+v.ddy_center())
    p.applyBC()
    a = 1/mesh.dx**2
    b = -1*(2./mesh.dx**2+2./mesh.dy**2)
    c = 1/mesh.dy**2
    w = para['w']
    max_iteration = para['max_teration_poission']
    tolerance = para['tolerance_poisson']
    status = False
    print "call poisson solver.\n"

    for n in range(max_iteration):
        residual = 0.
        for i in range(1,mesh.nx+1):
            for j in range(1,mesh.ny+1):
                p_old = p.value[j,i]
                #print "j,i:",j,i
                #print "p_old,",p_old
                #print "b:",b," a:",a," c:",c
                #print "Q(i,j):",  Q.value[j,i]
                #print "P(i+1,j):",p.value[j,i+1]
                #print "P(i-1,j):",p.value[j,i-1]
                #print "P(i,j+1):",p.value[j+1,i]
                #print "P(i,j-1):",p.value[j-1,i]
               
                p.value[j,i] = 1./b * (Q.value[j,i]-( a*(p.value[j,i+1]+p.value[j,i-1])+c*(p.value[j+1,i]+p.value[j-1,i])))
                #print "p(j,i): ",p.value[j,i]
                #print "p_old",p_old
                #print "p_new",p.value[j,i]
                p.value[j,i] = w * p.value[j,i] + (1-w) * p_old
                residual = max(residual, abs(p.value[j,i]-p_old))
        #print "iteration: ",n+1
        #print "    residual: ",residual
        if (residual < tolerance):
            print "poisson solver converged in ",n+1,"steps.\n"
            status = True
            p.applyBC()
            return status,p

    print "Poisson solver did not converge, exiting!"
    return status,p
    sys.exit()


            

