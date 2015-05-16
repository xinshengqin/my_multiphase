#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
import numpy as np

class mesh_1d(object):
    def __init__(self, xmin, xmax, delta1, deltan, name = 'default_mesh_1d'):
    #The four input arguments are coordinates of start point, coordinates of end point,
    #length of first cell and length of last cell
    #required number of cells are automatically computed
    #for uniform mesh (delta1==deltan), if (xmax-xmin)/delta1 != int
    #then actual delta1 = delta1, while actual deltan != deltan
    #The rule is delta1 is always guaranteed while deltan is not
        self.startPoint = float(xmin)
        self.endPoint = float(xmax)
        self.name = name
        self.length=xmax-xmin
        self.numberOfCells = self.get_n(delta1,deltan,self.length) #not include 4 ghost cells
        self.delta1=delta1
        self.deltan=deltan
        self.min_delta = min(delta1,deltan)
        #ratio between two adjacent cells, r = delta_{i+1} / delta_{i}  
        self.r=self.get_r_given_delta1(delta1,self.length,self.numberOfCells)
        self.edge=np.zeros(self.numberOfCells+1+4)#2 ghost cells at left and 2 at right
        self.cell_length = np.zeros(self.numberOfCells+4)
        for i in range(self.numberOfCells+4):
            self.cell_length[i] = delta1*self.r**(i-2)
        self.edge[2] = xmin
        self.edge[1] = self.edge[2] - self.cell_length[1]
        self.edge[0] = self.edge[1] - self.cell_length[0]
        for i in range(self.edge.size-3):
            self.edge[i+3] = self.edge[i+2]+self.cell_length[i+2]
        self.center = np.zeros(self.edge.size)
        self.center[1:] = (self.edge[:-1]+self.edge[1:])/2
        self.center[0] = self.edge[0] - self.cell_length[0]/self.r/2
        #print "self.edge: ",self.edge
        #print "self.length: ",self.cell_length
        #print "self.center: ", self.center


    def get_n(self,delta1,deltan,length,debug=0):
        #given delta1 (size of 1st cell), deltan (size of last cell) and length of the domain
        #giving it a guess of n can accelerate the computation
        #return number of cells, n
        #deltan may be slightly different from given value
        if debug > 0:
            print "call get_n()"

        deltan = float(deltan)
        delta1 = float(delta1)
        length = float(length)
        n = length*2/(delta1+deltan)#initial guess of n
        r1 = deltan/delta1
        if abs(deltan-delta1) < 1e-6:#deltan=delta1
            return int(round(n))#round n into a integer
        max_it = 100000
        convergence_rule = 0.0000001

        #newton-raphson iteration
        #lhs = delta1*(1-r1**(n/(n-1)))/(1-r1**(1/(n-1)))
        #solve for f =  a*(1-b^(x/(x-1)))/(1-b^(1/(x-1)))-L = 0
        f =  delta1*(1-r1**(n/(n-1))) / (1-r1**(1/(n-1)))- length
        from math import log
        for i in range(max_it):
            if debug > 0:
                print "iteration: ",i
                print "n:",n
                print "residual:",f
                print "\n"
            if abs(f) < convergence_rule:
                break
            #f_prime got from wolfram alpha
            f_prime = ( (delta1*  log(r1) )*(r1**(n/(n-1))-r1**(1/(n-1))) ) / ((n-1)**2 * (r1**(1/(n-1))-1)**2)
            n = n - f/f_prime #new n
            f =  delta1*(1-r1**(n/(n-1)))/(1-r1**(1/(n-1)))- length # new f

        if debug > 0 :
            print "total number of iteration for getting n:",i
            print "residual:",abs(f)
            print "n:",n

        return int(round(n))#round n into a integer

    def get_r_given_delta1(self,delta1,length,n,debug=0):
        #return ratio between two adjacent cells, r = delta_{i+1} / delta_{i}
        #solve delta1*(1-r^n)/(1-r) = length to get r
        if debug > 0 :
            print "\n"*3
            print "Given delta1, length of edge, number of cells along this edge"
            print "call get_r_given_delta1"

        delta1 = float(delta1)
        length = float(length)
        n = float(n)
        max_it = 2000
        convergence_rule = 0.0000001
        r = (length/n/delta1)**(1/(n-1))#initial value for r
        if abs(r - 1) < 1e-6:#uniform grid
            if debug > 0 :
                print "Warning! The inputs indicates a uniform grid!"
            i=0
        else:#nonuniform grid
            #newton-raphson iteration
            #solve for f(r) = (1-r^n)/(1-r) - L/delta1 = 0
            f = (1-r**n)/(1-r) - length/delta1 #f0
            for i in range(max_it):
                if debug == 2:
                    print "iteration:",i
                    print "r=",r
                    print "residual:",f
                if abs(f) < convergence_rule:
                    break
                f_prime = ((n-1)*r**(n+1)-n*r**n+r)/((1-r)**2*r)
                r = r - f/f_prime #new r
                f = (1-r**n)/(1-r) - length/delta1 #new f
            if i == max_it-1:
                print "Warning!! Maximum iteration reached. Something was wrong!"
        clusterRatio=r**(n-1)#ratio of size of last cell to first cell, delta_n/delta_1 

        if debug > 0:
            print "#"*80
            print "result summary"
            print "iteration times:", i
            print "ratio between two adjacent cells, r = delta_{i+1} / delta_{i} = ", r
            print "ratio of size of last cell to first cell, delta_n/delta_1 =",clusterRatio
            print "total number of cells along this edge, n =",n
            print "size of first cell:", delta1
            print "size of last cell:", delta1*r**(n-1)

        return r
