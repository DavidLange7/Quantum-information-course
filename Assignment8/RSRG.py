#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 18:12:08 2023

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
plt.style.use('ggplot')
import pylab as plb
import datetime

plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (18,12)

def tenspr(A,B):
    a1 = A.shape[0]
    a2 = A.shape[1]
    b1 = B.shape[0]
    b2 = B.shape[1]
    C = np.zeros((a1*b1,a2*b2))
    for m in range(0,C.shape[0]):
        for n in range(0,C.shape[1]):
            C[m,n] = A[m//b1,n//b2]*B[m%b1,n%b2]
    return C

def stp1(sides, lambd):
    
    pau_x = np.array([[0,1], [1,0]])
    pau_z = np.array([[1,0], [0,-1]])


    H = 0
    
    for i in range(1,sides+1):
        if i==1:
            tmp = pau_z
        else: 
            tmp = np.identity(2)
            
        for j in range(2,sides+1):
            if (j==i and i != 1):
                tmp = tenspr(tmp, pau_z)
            else:
                tmp = tenspr(tmp, np.identity(2))
        H = H + tmp
    H = lambd*H
    
    H2 = 0
    
    for i in range(1,sides):
        bool1 = False
        
        if i==1:
            tmp2 = pau_x
        else:
            tmp2 = np.identity(2)
            
        for j in range(2,sides+1):
            if j==i and i!=1:
                tmp2 = tenspr(tmp2, pau_x)
                bool1 = True
            elif bool1==True:
                tmp2 = tenspr(tmp2, pau_x)
                bool1 = False
            elif i==1 and j==2:
                tmp2 = tenspr(tmp2, pau_x)
                bool1 = False
            else:
                tmp = tenspr(tmp2, np.identity(2))
        H2 = H2 + tmp2
        
    H = H + H2
        
    H_int_L = np.identity(2)
    H_int_R = pau_x
        
    for i in range(2,sides+1):
        if i==sides:
            H_int_L = tenspr(H_int_L, pau_x)
        else:
            H_int_L = tenspr(H_int_L, np.identity(2))
        H_int_R = tenspr(H_int_R, np.identity(2))
                
    H_int = tenspr(H_int_L, H_int_R)
    
    return(H_int, H_int_L, H_int_R, H)

def H_stp(H, sides):
    
    H2R = np.identity(2)
    
    for i in range(1,sides+1):
        if i==1:
            H2L = tenspr(H, np.identity(2))
        else: 
            H2L = tenspr(H2L, np.identity(2))
        
        if i==sides:
            H2R = tenspr(H2R, H)
        else: 
            H2R = tenspr(H2R, np.identity(2))
            
    H_new = H2L + H2R
    
    return(H_new)

def interaction(H_int_L, H_int_R, sides):
    
    tmp = np.identity(2)
    
    for i in range(1,sides+1):
        if i==sides:
            tmp = tenspr(tmp, H_int_L)
        else: 
            tmp = tenspr(tmp, np.identity(2))
        
        H_int_R = tenspr(H_int_R, np.identity(2))
        
    H_int_L = tmp
    
    return(H_int_L, H_int_R)

def RSRG(lambd, sides, M):
    #print(datetime.datetime.now())
    H_int, H_int_L, H_int_R, H = stp1(sides, lambd)
        
    n = sides
    gs_e = []
    
    
    for i in range(1,M+1):
        n = 2*n
        print(n)

        H_new = H_stp(H, sides)
        #print(np.shape(H_int))
        
        H = H_new + H_int
        w, v = np.linalg.eigh(H)
        l = np.argsort(w)
        
        w = w[l]
        v = v[l]


        H_int_L, H_int_R = interaction(H_int_L, H_int_R, sides)

        
        P = v[:, 0:2**sides]
        P_d = P.conj().T
        
        H = np.matmul(P_d, np.matmul(H, P))
        H_int_L = np.matmul(P_d, np.matmul(H_int_L,P))
        H_int_R = np.matmul(P_d, np.matmul(H_int_R,P))
        
        H_int = tenspr(H_int_L, H_int_R)

        gs_e.append(w[0]/n)
    
    #plt.figure(1)
    #plt.ylabel('E/N')
    #plt.xlabel('iteration step')

    #plt.plot(gs_e, 'o:r', label = 'RSRG')
    #plt.savefig('gs_iter.pdf', dpi=1000,bbox_inches='tight')
    
    #print(datetime.datetime.now())
    return(gs_e)
    

def ch_lambd(sides, M, maxlambd):
    
    gs_res = []
    
    lambd_space = np.linspace(0,maxlambd,100)
    for i in range(len(lambd_space)):
        tmp = RSRG(lambd_space[i], sides, M)
        gs_res.append(tmp[-1])
    '''
    plt.figure(1)
    plt.ylabel('E/N')
    plt.xlabel('\u03BB')

    plt.plot(lambd_space, gs_res, label = 'RSRG')

    spline = UnivariateSpline(lambd_space, gs_res, s = 0, k = 4)
    
    mf_lmdspc1 = lambd_space[0:np.where(lambd_space==2)[0][0]]
    mf_lmdspc2 = lambd_space[np.where(lambd_space==2)[0][0]-1:-1]

    plt.plot(mf_lmdspc1, -1-mf_lmdspc1**2/4, color = 'blue', label = 'MF')
    plt.plot(mf_lmdspc2, - np.abs(mf_lmdspc2), color = 'blue')

    d2spl = spline.derivative(n=2)
    plt.legend()
    
    plt.figure(2)
    plt.ylabel('d²(E/N)/d\u03BB²')
    plt.xlabel('\u03BB')
    plt.plot(lambd_space,d2spl(lambd_space))
    
    fig, axs = plt.subplots(2)

    axs[0].plot(lambd_space, gs_res, 'ro', label = 'RSRG')
    axs[0].plot(lambd_space, spline(lambd_space), color = 'black', alpha = 0.7, label = 'spline')
    axs[0].set_ylabel('E/N')
    #axs[0].xlabel('\u03BB')
    axs[0].plot(mf_lmdspc1, -1-mf_lmdspc1**2/4, color = 'blue', label = 'MF')
    axs[0].plot(mf_lmdspc2, - np.abs(mf_lmdspc2), color = 'blue')
    axs[0].legend()
    axs[1].plot(lambd_space,d2spl(lambd_space))
    axs[1].set_ylabel('d²(E/N)/d\u03BB²')
    axs[1].set_xlabel('\u03BB')
    #plt.savefig('RSRG_MF.pdf', dpi=1000,bbox_inches='tight')
    '''
