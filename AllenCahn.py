#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 17:04:48 2016
Allen-Cahn style demixing
@author: dyanni3
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
D = 0.5
size = 300  # size of the 2D grid
dx = 2./size  # space step
T = 0.5  # total time
dt = .9 * dx**2/2  # time step
n = int(T/dt) #number of time slices
U = 2*np.random.rand(size, size)-1 #The domain
#X,Y=np.ogrid[-int(size/2):int(size/2),-int(size/2):int(size/2)]
#U=np.ones((size,size))*(X**2+Y**2<=400)
#U=2*U-1

def laplacian(Z):
    Ztop = Z[0:-2,1:-1]
    Zleft = Z[1:-1,0:-2]
    Zbottom = Z[2:,1:-1]
    Zright = Z[1:-1,2:]
    Zcenter = Z[1:-1,1:-1]
    return (1.1*Ztop + Zleft + Zbottom + Zright - 4 * Zcenter) / dx**2
#%%
im1=plt.imshow(U, cmap=plt.cm.bwr, extent=[-1,1,-1,1]);
#%%
def update(U,t):
    # We simulate the PDE with the finite difference method.
    for i in range(int(t*n/10000)):
        # We compute the Laplacian of u and v.
        deltaU = laplacian(U)
        # We take the values of u and v inside the grid.
        Uc = U[1:-1,1:-1]
        # We update the variables.
        U[1:-1,1:-1] = Uc + dt * (D * deltaU + -(Uc**3)+(2*Uc)) 
        # Neumann conditions: derivatives at the edges
        # are null.
        for Z in (U,U):
            Z[0,:] = Z[1,:]
            Z[-1,:] = Z[-2,:]
            Z[:,0] = Z[:,1]
            Z[:,-1] = Z[:,-2]
    return U
#%%    
im2=plt.imshow(U, cmap=plt.cm.seismic, extent=[-1,1,-1,1]);

#%%
def drop(U,radius,intensity):
    unpadded_size=U.shape[0]-2
    u=np.random.randint(0,(U.shape[0]-2)**2)
    x0=int(u/unpadded_size)+1
    y0=(u%unpadded_size)+1
    for i in range(U.shape[0]):
        for j in range(U.shape[1]):
            if (((i-x0)**2) + ((j-y0)**2)<=radius**2):
                U[i,j]=intensity
    return 
#%%
V=U.copy()
for i in range(len(U)):
    for j in range(len(U)):
        if U[i,j]<=0:
            V[i,j]=-1
        else:
            V[i,j]=1

im3=plt.imshow(V, cmap=plt.cm.bwr, extent=[-1,1,-1,1]);           
#%%
from PIL import Image
def data2RGB(x):
	if x==0.0:
		return(0,0,0);
	elif x>0:
		return(255,0,0)
	else:
		return(0,0,255)
l=list(map(data2RGB,np.ravel(U)))
#img_layer=np.array(l).reshape(self.lattice.shape[1],self.lattice.shape[1],3)
im = Image.new('RGB', [U.shape[1],U.shape[0]]);
im.putdata(l);
im.show()

#%%
def findTime(radius,strength):
    D = 0.25
    size = 200  # size of the 2D grid
    dx = 2./size  # space step
    T = 0.5  # total time
    dt = .9 * dx**2/2  # time step
    n = int(T/dt) #number of time slices
    #U = 2*np.random.rand(size, size)-1 #The domain
    X,Y=np.ogrid[-int(size/2):int(size/2),-int(size/2):int(size/2)]
    U=np.ones((size,size))*(X**2+Y**2<=radius**2)
    U=2*U-1
    t=0
    while(U[U>0].size>0):
        t+=1
        for i in range(int(n/1000)):
            # We compute the Laplacian of u and v.
            deltaU = laplacian(U)
            # We take the values of u and v inside the grid.
            Uc = U[1:-1,1:-1]
            # We update the variables.
            U[1:-1,1:-1] = Uc + dt * (D * deltaU + -(Uc**3)+(2*Uc)-strength) 
            # Neumann conditions: derivatives at the edges
            # are null.
            for Z in (U,U):
                Z[0,:] = Z[1,:]
                Z[-1,:] = Z[-2,:]
                Z[:,0] = Z[:,1]
                Z[:,-1] = Z[:,-2]
    return(t)
#%%
def findRadius(U):
    h=int(U.shape[0]/2)
    return len(np.where(U[h,:]>0)[0])
    
def rVt(U,strength):
    r=[]
    while(U[U>0].size>0):
        for i in range(int(n/1000)):
            # We compute the Laplacian of u and v.
            deltaU = laplacian(U)
            # We take the values of u and v inside the grid.
            Uc = U[1:-1,1:-1]
            # We update the variables.
            U[1:-1,1:-1] = Uc + dt * (D * deltaU + -(Uc**3)+(2*Uc)-strength) 
            # Neumann conditions: derivatives at the edges
            # are null.
            for Z in (U,U):
                Z[0,:] = Z[1,:]
                Z[-1,:] = Z[-2,:]
                Z[:,0] = Z[:,1]
                Z[:,-1] = Z[:,-2]
            r.append(findRadius(U))
    return r
    
#%%
def video(lattice,tlist):
    for (x,t) in enumerate(tlist):
        if x<10:
            plt.imsave('XartVid00%d.png'%x,lattice,cmap=plt.cm.Blues)
        elif x<100:
            plt.imsave('XartVid0%d.png'%x,lattice,cmap=plt.cm.Blues)
        else:
            plt.imsave('XartVid%d.png'%x,lattice,cmap=plt.cm.Blues)
        lattice=update(lattice,t)
        if (np.random.rand()>.5 and x>50):
            drop(lattice,t*np.random.rand()**2,np.mean(U))
    return