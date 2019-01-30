# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
import numpy as np
from scipy.special import kn #Bessel
from matplotlib import pyplot as plt

#%%
def beta(x,y,D,gamma):
    r2 = x**2 + y**2
    num = np.sqrt(D*r2*(4*(D**2)*gamma))
    denom = 4*D**2
    return(num/denom)

def q(x,y,alpha,D,gamma):
    pref = alpha/(2*np.pi*D)
    bessel = kn(0,2*beta(x,y,D,gamma))
    return(pref*bessel)
    
#%% # set up the domain
x = np.linspace(-1,1,100)
y = np.linspace(-1,1,100)
xx, yy = np.meshgrid(x,y,sparse=True)
plt.imshow(5*q(xx,yy,1,1,1)+3*q(xx-.1,yy-.1,1,1,1)+30*q(xx+1.1,yy-.1,1,1,1))

#%% #make a random initial rho
from scipy import signal
cell_initial_pos = tuple([tuple(np.random.randint(0,100,2)) for i in range(300)])
rho_initial = np.stack([signal.unit_impulse((100,100),i) for i in cell_initial_pos])
rho_initial = np.sum(rho_initial,axis=0)
plt.imshow(rho_initial,cmap=plt.cm.gray)

#%%
qtot = np.zeros((100,100))
for i in range(100):
    for j in range(100):
        if j%90==0 and i%10==0:
            print(i)
        this_q = q(xx-x[j],yy-y[i],10,.1,20)
        this_q[i,j]=q(.01,.01,10,.1,20)
        qtot += rho_initial[i,j]*this_q
plt.imshow(qtot)

#%%
atot = np.zeros((100,100))
rho_a = rho_initial*(qtot>80)
for i in range(100):
    for j in range(100):
        if j%90==0 and i%10==0:
            print(i)
        this_a = q(xx-x[j],yy-y[i],100,1,.1)
        this_a[i,j] = q(.01,.01,100,1,.1)
        atot += rho_a[i,j]*this_a
plt.imshow(atot,cmap=plt.cm.jet)

#%%
rho_safe = rho_initial*(atot>700)
plt.imshow(rho_safe)
#%%
rho_die = rho_initial*(atot<=700)
D_kernel = np.array([[0,1,0],[1,-4,1],[0,1,0]])
rho_initial = rho_initial + .1*signal.convolve2d(rho_initial,D_kernel,mode='same')
rho_initial = rho_initial + rho_safe*.1*(10-rho_initial)/10
rho_initial = rho_initial - .5*rho_die
plt.imshow(rho_initial,cmap=plt.cm.gray)

#%%
import sys
def progressBar(value, endvalue, bar_length=20):
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

def full_sim(alpha_q=10, D_q=.1, gamma_q=20,
             alpha_a=100, D_a=1, gamma_a=20,
             q_thresh=80, a_thresh=700,
             D_rho=.05, dt=.05, K=10, rgrow=10, rkill=.5,
             rho_initial='none',nsteps=1):
       for n in range(nsteps):
           print(" ")
           print("Simulation running step %d out of %d"%(n+1,nsteps))
           if(rho_initial=='none'):
               cell_initial_pos = tuple([tuple(np.random.randint(0,100,2)) for i in range(300)])
               rho_initial = np.stack([signal.unit_impulse((100,100),i) for i in cell_initial_pos])
               rho_initial = np.sum(rho_initial,axis=0)
           rho_fig = plt.figure()
           plt.imshow(rho_initial,cmap=plt.cm.gray,vmin=0,vmax=10)
           path = "C:\\Users\\dyanni3\\Desktop\\santola_sims\\RHO\\prhofig%d.png"%(n-1)
           rho_fig.savefig(path)
           qtot = np.zeros((100,100))
           print("Configuring QS concentrations")
           for i in range(100):
               for j in range(100):
                   progressBar((100*i)+j,10000)
                   this_q = q(xx-x[j],yy-y[i],alpha_q,D_q,gamma_q)
                   this_q[i,j]=q(.01,.01,alpha_q,D_q,gamma_q)
                   qtot += rho_initial[i,j]*this_q
           qs_fig = plt.figure()
           plt.imshow(qtot)
           path = "C:\\Users\\dyanni3\\Desktop\\santola_sims\\QS\\pqsfig%d.png"%n
           qs_fig.savefig(path)
           atot = np.zeros((100,100))
           rho_a = rho_initial*(qtot>q_thresh)
           print(" ")
           print("Configuring Antibiotic concentrations")
           for i in range(100):
               for j in range(100):
                   progressBar((100*i)+j,10000)
                   this_a = q(xx-x[j],yy-y[i],alpha_a,D_a,gamma_a)
                   this_a[i,j] = q(.01,.01,alpha_a,D_a,gamma_a)
                   atot += rho_a[i,j]*this_a
           a_fig = plt.figure()
           plt.imshow(atot,cmap=plt.cm.jet)
           path = "C:\\Users\\dyanni3\\Desktop\\santola_sims\\AB\\pabfig%d.png"%n
           a_fig.savefig(path)
           rho_safe = rho_initial*(atot>a_thresh)
           plt.imshow(rho_safe)
           rho_die = rho_initial*(atot<=a_thresh)
           D_kernel = np.array([[0,1,0],[1,-4,1],[0,1,0]])
           rho_initial = rho_initial + dt*D_rho*signal.convolve2d(rho_initial,D_kernel,mode='same')
           rho_initial = rho_initial + dt*rho_safe*rgrow*(K-rho_initial)/K
           rho_initial = rho_initial - dt*rkill*rho_die
           rho_fig = plt.figure()
           plt.imshow(rho_initial,cmap=plt.cm.gray,vmin=0,vmax=10)
           path = "C:\\Users\\dyanni3\\Desktop\\santola_sims\\RHO\\prhofig%d.png"%n
           rho_fig.savefig(path)
       return(rho_initial)