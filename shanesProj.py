# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 14:57:31 2019

@author: dyanni3
"""

#%%
import numpy as np
from scipy.linalg import circulant
from scipy.optimize import minimize
from matplotlib import pyplot as plt

#%%
def build_c(beta, N, mode = 'bipartite'):
    
    if mode=='bipartite':
        A = circulant([i%2 for i in range(1,N+1)])
        one = np.ones(N)
        t1 = beta*(A*np.outer(1/np.dot(one,A),one))
        t2 = (1-beta)*np.eye(N)
        return(t1+t2)
    if mode=='neighbor':
        l = np.zeros(N)
        l[:3] = 1
        A = circulant(l)
        one = np.ones(N)
        t1 = beta*(A*np.outer(1/np.dot(one,A),one))
        t2 = (1-beta)*np.eye(N)
        return(t1+t2)
    if mode=='full':
        A = np.ones((N,N))
        one = np.ones(N)
        t1 = beta*(A*np.outer(1/np.dot(one,A),one))
        t2 = (1-beta)*np.eye(N)
        return(t1+t2)
        
def W(v,c, alpha):
    right = (1-v)**alpha
    right = np.dot(c,right)
    return(-np.dot(v**alpha, right))
    
def grad_W(v,c, alpha):
    left1 = alpha*(v**(alpha-1))
    right1 = np.dot(c,(1-v)**alpha)
    t1 = left1*right1
    left2 = alpha*(1-v)**(alpha-1)
    right2 = np.dot(c,v**alpha)
    t2 = left2*right2
    return(t1-t2)
    
def W_and_grad(v, c=build_c(1,10),alpha=1):
    return( W(c,alpha,v), grad_W(c,alpha,v))
    
#%%
def spec_optimized(beta=1, N=10, mode='bipartite', alpha=1, max_iter = 100):
    bnds = tuple([(0,1) for i in range(N)])
    c = build_c(beta,N,mode)
    best = 999
    for i in range(max_iter):
        res = minimize(W, np.random.random(N), method='nelder-mead', jac=grad_W,
                   options={'disp': False, 'eps':1e-12, 'gtol':1e-9, 'ftol':1e-9},
                   args = (c,alpha), bounds=bnds)
        if res.fun<best:
            best = res.fun
            v = res.x
    else:
        res = minimize(W, np.random.random(N), method='L-BFGS-B', jac=grad_W,
                   options={'disp': True}, args = (c,alpha), bounds=bnds)
        pass
    spec = np.average([2*(max([vi, (1-vi)])-.5) for vi in v])
    #spec = np.average(v)
    return(spec,best)
    
#%% parameter sweep
def make_plot(mode='bipartite'):
    alphas = np.linspace(0.01,1.5,50)
    betas = np.linspace(0,1,4)
    specs = [[] for i in range(10)]
    for i,beta in enumerate(betas):
        print(i)
        for j,alpha in enumerate(alphas):
            print(j)
            s,res = spec_optimized(beta=beta, alpha=alpha, mode = mode ,N=10)
            specs[i].append(s)
            
    fig = plt.figure(figsize = (10,10))
    plt.xlabel("Specialization Power",fontsize = 14)
    plt.ylabel("Specialization",fontsize=14)
    plt.title("%s network"%mode,fontsize=14)
    for i in range(4):
        plt.plot(alphas,specs[i],marker='o',lw=0, label = betas[i])
    plt.legend()
    return(fig)
#%%
def gradient_descent(eps, v0, beta=1, N=10, mode='bipartite', alpha=1, max_iter = 1000):
    ws = []
    c = build_c(beta,N,mode)
    for i in range(max_iter):
        v0 = v0 + eps*grad_W(v0,c,alpha)
        v0[v0>1]=1
        v0[v0<0]=0
        ws.append(W(v0,c,alpha))
    return(v0,ws)
    
def spec(v):
    return(np.average([2*(max([vi, (1-vi)])-.5) for vi in v]))