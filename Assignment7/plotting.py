#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 17:03:03 2022

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pylab as plb

#plb.rcParams['font.size'] = 55
#plt.rcParams["figure.figsize"] = (22,14)

plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (18,12)

#%%

E = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/eigenvalues.txt')#, delimiter = ',')
psi = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/eigenvectors.txt')


l = int(np.sqrt(len(psi)))
psi = np.split(psi, l)
psi2 = [psi[i]*psi[i] for i in range(l)]
    
#%% THIS IS FOR THE FIRST 4 PLOTS TO SEE FIRST K = 128 LEVELS (ITS A N=7 SYSTEM). LAMBDA HAS 5 DIFFERENT VALUES HERE
E_ar = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/eigenvalues_deg.txt')#, delimiter = ',')

lbls = ['\u03BB = 0', '\u03BB = 0.1', '\u03BB = 1', '\u03BB = 3', '\u03BB = 10']
for i in range(1,5):
    plt.figure(i)
    plt.plot(E_ar[i], '.', label = lbls[i],marker = 'o', markersize = 9)
    #plt.ylabel('E')
    plt.xlabel('Level k')
    plt.legend()
    lbls1 = ['L0.pdf', 'L0_1.pdf', '1.pdf', '3.pdf', '10.pdf']
    plt.savefig(lbls1[i], dpi=1000,bbox_inches='tight')
    
#%% THIS IS FOR CHECKING THE FIRST FEW LEVELS AS A FUNCTION OF LAMBDA SO LAMBDA GOES FROM 0-10 in 300 steps
E_arr = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/eigenvalues_gsdeg.txt')#, delimiter = ',')

E_arr = np.transpose(E_arr)
lmspc = [(i*1.0)/30 for i in range(0,300)]

#plt.xlim([0,10])
'''
plt.plot(lmspc, E_arr[0])
plt.plot(lmspc, E_arr[1])
plt.plot(lmspc, E_arr[2])
plt.plot(lmspc, E_arr[4])
'''
plt.figure()
plt.ylabel('E')
plt.xlabel('\u03BB')
for i in range(0,len(E_arr)):
    plt.plot(lmspc, E_arr[i], color = 'gray', linewidth = 0.8)
plt.savefig('deglambda.pdf', dpi=1000,bbox_inches='tight')
#%%sparsity stuff
sps = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/sparsity.txt')#, delimiter = ',')
res = []
plt.figure('sparse')
plt.xlabel('Number of sides N')
plt.ylabel('sparsity')
for i in range(1,13):
    res.append(((2**i*2**i)-sps[i-1])/(2**i*2**i))
    plt.plot(i,res[i-1], '*', color = 'black', markersize = 12)
plt.savefig('sparcity.pdf', dpi=1000,bbox_inches='tight')
#%% this is for the magnetization stuff
M = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/magnetization.txt')#, delimiter = ',')
gs = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment7/gs.txt')#, delimiter = ',')

res = np.matmul(gs,M)
res = np.matmul(res,gs)

#sides =7 for this one..
M = np.array([0.0014070556555936115,0.023537022207499737, 0.056709500846935226, 0.37488205687424114, 0.5751460457680815, 0.7247498489956841, 0.8967415155816775, 0.9644389095511604, 0.9971312439057464, 0.9999709718722627, 0.9999997093876648])
lmd= np.array([0.0005, 0.01, 0.1, 0.5, 0.75, 1, 1.5, 3, 10, 100, 1000])
fig, ax = plt.subplots()
ax.plot(lmd, M)
ax.plot(lmd, M, 'bo', label = 'N=7')
plt.xlabel('\u03BB')
plt.ylabel('M')
ax.set_xscale('log')
plt.savefig('magnet.pdf', dpi=1000,bbox_inches='tight')
