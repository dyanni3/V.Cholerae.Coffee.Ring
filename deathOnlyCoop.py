#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 10:51:05 2016
2D IBM
@author: dyanni3
"""

import numpy as np
class Lattice(object):
    mask=np.array([[0,1,0],[1,0,1],[0,1,0]])
    colors=np.array(['r','b'])
    def __init__(self,size=5,radius=0):
        x,y=size,size;
        if not radius:
            self.lattice=np.random.choice(['r','b','0'],size=(x,y),p=[(.475),(.475),(.05)])
        else:
            X,Y=np.ogrid[-int(size/2):int(size/2),-int(size/2):int(size/2)]
            U=np.ones((size,size))*(X**2+Y**2<=radius**2)
            U=2*U-1
            def U2rgb(num):
                num+=1
                return ['r','0','b'][int(num)]
            self.lattice=np.array(list(map(U2rgb,np.ravel(U)))).reshape(x,y)
            
        return
        
    def evolve(self,n_steps):
        for t in range(n_steps):
            
            #pick a site (i,j) 
            unpadded_size=self.lattice.shape[0]-4
            u=np.random.randint(0,(self.lattice.shape[0]-4)**2)
            i=int(u/unpadded_size)+2
            j=(u%unpadded_size)+2

            #24 member neighborhood (for cooperative benefits counting)
            neighb5x=np.ravel(self.lattice[i-2:i+3,j-2:j+3])
            #4 member neighborhood (for doubling into empty sites)
            neighb3x=self.lattice[i-1:i+2,j-1:j+2]
  
            altruist_count=len(np.where(neighb5x=='b')[0])
            growth_bonus_blue=(altruist_count**10)/(12**10)
            growth_bonus_red=10
            
            #respond to neighborhood
            if self.lattice[i,j]=='r':
                for (row,entries) in enumerate(neighb3x):
                    for (column,item) in enumerate(entries):
                        if item=='0':
                            if growth_bonus_red*growth_bonus_blue*np.random.random()>.85:
                                self.lattice[i+row-1,j+column-1]='r'
                if np.random.random()>.95:
                    self.lattice[i,j]='0'
            elif self.lattice[i,j]=='b':
                for (row,entries) in enumerate(neighb3x):
                    for (column,item) in enumerate(entries):
                        if item=='0':
                            if growth_bonus_blue*np.random.random()>.85:
                                self.lattice[i+row-1,j+column-1]='b'
                if np.random.random()>.95:
                    self.lattice[i,j]='0'
                    
        return
    def view(self):
        from PIL import Image
        def data2RGB(x):
            if x=='0':
                return(0,0,0);
            elif x=='r':
                return(255,0,0)
            elif x=='g':
                return(0,255,0)
            elif x=='b':
                return(0,0,255)
        l=list(map(data2RGB,np.ravel(self.lattice[:,:])))
        im = Image.new('RGB', [self.lattice.shape[1],self.lattice.shape[0]]);
        im.putdata(l);
        im.show
        return im 
  #%%
def video(lattice):
    for t in range(50):
        im=lattice.view()
        im.save('coopVid10x%d.png'%t)
        lattice.evolve(100000)
    return
    
#%%
def plt_blueVred(lattice):
    nr=[]
    nb=[]
    for t in range(50):
        nr.append(len(np.where(lattice.lattice=='r')[0]))
        nb.append(len(np.where(lattice.lattice=='b')[0]))
        print(t)
        lattice.evolve(100000)
    return nr,nb
        
        