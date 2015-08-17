#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


filename = 'final'
path='6.0'
Nx=40
Ny=40

Lx=1.
Ly=1.
dx = Lx/Nx
dy = Ly/Ny
#U_inf = 1
#nu=0.01
#gamma = 0.66
#L=4.
#u=np.zeros(Ny,Nx)
#v=np.zeros(Ny,Nx)
#p=np.zeros(Ny,Nx)
u=[]
with open(path+'/c.csv','r') as inputfile:
    for line in inputfile:
        if line.startswith('#'):
            continue
        line=line.strip(' ').strip('\n').split(',')
        u.append(line[1])

u=np.array(u)
u=u.astype(np.float)


c_exact=[]
with open('0.0/c.csv','r') as inputfile:
    for line in inputfile:
        if line.startswith('#'):
            continue
        line=line.strip(' ').strip('\n').split(',')
        c_exact.append(line[1])

c_exact=np.array(c_exact)
c_exact=c_exact.astype(np.float)

error = 0
for i in range(u.size):
    error = error+dx*dy*abs(u[i]-c_exact[i])
print 'error: ',error

delta = u-c_exact
delta=delta.reshape(Ny+2,Nx+2)#include ghost cells

x=np.linspace(0,Lx+dx,Nx+2)
y=np.linspace(0,Ly+dy,Ny+2)

#zeta_edge=np.linspace(0,Ly,Ny+1)
#zeta_edge=zeta_edge[1:]
#zeta_center=zeta_edge-Ly/Ny/2
#y_edge = Ly * (1.0 - np.tanh(gamma * (Ly - zeta_edge)) / np.tanh(gamma*Ly))
#y_center = Ly * (1.0 - np.tanh(gamma * (Ly - zeta_center)) / np.tanh(gamma*Ly))

X,Y=np.meshgrid(x,y)
        
        




#fig, axs = plt.figure()
fig = plt.figure()
ax = fig.add_subplot(111)
#levels = np.linspace(-0.1, 1.2, 40)
#levels = np.linspace(-0.1, 1.2, 40)
#cs = plt.contourf(X, Y, p, levels=levels)
cs = plt.contourf(X, Y, delta,50,cmap=plt.cm.get_cmap('jet'))#50 specifys how many levels of colors is there in the plot
#cs2 = plt.contour(cs, levels=cs.levels[::4], colors = 'k', origin='upper', hold='on')
#fig.colorbar(cs, ax=axs[0], format="%.2f")
plt.colorbar(cs,format="%.2f",orientation='horizontal')
legend = plt.legend(loc='upper center', shadow=True,fontsize='x-large')
#plot mesh grid
for j in range(y.size):
    yg = np.zeros(x.size)+y[j]
    plt.plot(x,yg,'b-',linewidth=0.1)
for i in range(x.size):
    xg = np.zeros(y.size)+x[i]
    plt.plot(xg,y,'b-',linewidth=0.1)
plt.xlabel(r'x')
plt.ylabel(r'y')
ax.set_aspect('equal')
if not('myplot' in os.listdir('./')):
    os.mkdir('myplot')
fig.savefig('./myplot/contour_errorC_'+filename+'.png',bbox_inches='tight',dpi=200)

#cs = axs[1].contourf(X, Y, zdata, levels=[-1,0,1])
#fig.colorbar(cs, ax=axs[1])

plt.show()





