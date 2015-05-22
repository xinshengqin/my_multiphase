#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
import numpy as np
import os

class Mesh_1d(object):
    def __init__(self, xmin, xmax, delta1, deltan, name = 'default_mesh_1d'):
        #The four input arguments are coordinates of start point, coordinates of end point,
        #length of first cell and length of last cell
        #required number of cells are automatically computed
        #for uniform mesh (delta1==deltan), if (xmax-xmin)/delta1 != int
        #then actual delta1 = delta1, while actual deltan != deltan
        #The rule is delta1 is always guaranteed while deltan is not
        self.startPoint = float(xmin)
        self.endPoint = float(xmax)
        self.name = str(name)
        self.length=xmax-xmin
        self.numberOfCells = self.get_n(delta1,deltan,self.length) #not include 4 ghost cells
        self.delta1=float(delta1)
        self.deltan=float(deltan)
        self.min_delta = min(self.delta1,self.deltan)
        #ratio between two adjacent cells, r = delta_{i+1} / delta_{i}  
        self.r=self.get_r_given_delta1(delta1,self.length,self.numberOfCells)
        #no. of ghost cells = 2
        self.edge=np.zeros(self.numberOfCells+2)
        self.cell_length = np.zeros(self.numberOfCells+2)
        for i in range(self.numberOfCells+2):
            self.cell_length[i] = delta1*self.r**(i-1)
        self.edge[0] = float(xmin)
        #self.edge[1] = self.edge[2] - self.cell_length[1]
        #self.edge[0] = self.edge[1] - self.cell_length[0]
        for i in range(1,self.edge.size):
            self.edge[i] = self.edge[i-1] + self.cell_length[i]
        self.center = np.zeros(self.edge.size)
        self.center[1:] = (self.edge[:-1]+self.edge[1:])/2
        self.center[0] = self.edge[0] - self.cell_length[0]/2
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


    def plot_edgeVsIndex(self):
        index=np.arange(0,self.numberOfCells+2,1)
        plt.figure()
        plt.plot(index,self.edge,'ro')
        ax = plt.gca()
        ax.grid(b=True, which='major', color='b', linestyle='-')
        ax.grid(b=True, which='minor', color='r', linestyle='--')
        plt.minorticks_on()
        plt.xlabel('index i')
        plt.ylabel('edge coordinate')
        if not ('myplot' in os.listdir('./')):
            os.mkdir('./myplot')
        plt.savefig('./myplot/plot_edgeVsIndex_'+self.name+'.png', bbox_inches='tight')
        plt.close()

    def plot_centerVsIndex(self):
        index=np.arange(0,self.numberOfCells+2,1)
        plt.figure()
        plt.plot(index,self.center,'ro')
        ax = plt.gca()
        ax.grid(b=True, which='major', color='b', linestyle='-')
        plt.minorticks_on()
        ax.grid(b=True, which='minor', color='r', linestyle='--')
        plt.xlabel('index i')
        plt.ylabel('center coordinate')
        if not ('myplot' in os.listdir('./')):
            os.mkdir('./myplot')
        plt.savefig('./myplot/plot_centerVsIndex_'+self.name+'.png', bbox_inches='tight')
        plt.close()


class Mesh_2d(object):
    def __init__(self,x,y,name='default_mesh_2d'):
        #x,y should be two mesh_1d objects
        vedge_raw = np.meshgrid(x.edge,y.center)
        self.vedge_x = vedge_raw[0]#x coordinates of vertical edges 
        self.vedge_y = vedge_raw[1]#y coordinates of vertical edges 

        hedge_raw = np.meshgrid(x.center,y.edge)
        self.hedge_x = hedge_raw[0] #x coordinates of horizontal edges
        self.hedge_y = hedge_raw[1] #y coordinates of horizontal edges

        center_raw =  np.meshgrid(x.center,y.center)
        self.center_x = center_raw[0]#x coordinates of centers of cells
        self.center_y = center_raw[1]#y coordinates of centers of cells
        
        self.dx = x.min_delta
        self.dy = y.min_delta
        self.min_delta = min(x.min_delta,y.min_delta)
        self.nx = x.numberOfCells#number of values in x direction, not including ghost cells
        self.ny = y.numberOfCells
        self.lx = x.length
        self.ly = y.length
        self.x0 = x.startPoint
        self.xn = x.endPoint
        self.y0 = y.startPoint
        self.yn = y.endPoint
        self.name = name



        #print "x coordinates of vertical edges: ",self.vedge_x
        #print "y coordinates of vertical edges: ",self.vedge_y
        #print "x coordinates of horizontal edges: ",self.hedge_x
        #print "y coordinates of horizontal edges: ",self.hedge_y
        #print "x coordinates of centers of cells: ",self.center_x
        #print "y coordinates of centers of cells: ",self.center_y
    def write(self):
        file = open(self.name+'.csv','w')
        file.write('#No.,vedge_x,vedge_y,hedge_x,hedge_y,center_x,center_y\n')
        vx=np.ravel(self.vedge_x)
        vy=np.ravel(self.vedge_y)
        hx=np.ravel(self.hedge_x)
        hy=np.ravel(self.hedge_y)
        cx=np.ravel(self.center_x)
        cy=np.ravel(self.center_y)
        xy=np.vstack((vx,vy,hx,hy,cx,cy))
        xy=np.transpose(xy)
        np.set_printoptions(formatter={'float': lambda x: format(x, '+9.6E')})
        for i,item in enumerate(xy):
            line = '{:>9s}'.format(str(i))+','+",".join(str(item).lstrip('[').rstrip(']').split())+'\n'
            #line = str(i)
            #for j,item2 in enumerate(item):
            #    line = line+',{:16s}'.format(str(item2)) 
            #line = line+'\n'
            file.write(line)
        

