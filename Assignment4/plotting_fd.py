#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:39:57 2022

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pylab as plb
plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (14,16)
#plt.rcParams["figure.figsize"] = (16,10)

#%%
Z = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment4/Z.txt')#, delimiter = ',')
E = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment4/E.txt')#, delimiter = ',')
xmesh = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment4/mesh.txt')
#%%
plt.figure('Eigenvalue errors_2')
for i in range(0,10):
    plt.plot(xmesh, Z[i]+E[i])
#plt.plot(xmesh, wf_r*(max(wf_l)/max(wf_r)), 'gx', alpha = 0.3)
plt.plot(xmesh, 0.5*xmesh**2, color = 'black',  linestyle ='dashed')
E_analytic = np.array([(n+(1/2)) for n in range(0,len(E))])
plt.xlabel('X')
plt.ylabel('Energy')
plt.ylim([0,10])
plt.fill_between(xmesh, 0.5*xmesh**2,color='g',alpha=.4)
#plt.savefig('assgn4_wfs_elong.pdf', dpi=1000,bbox_inches='tight')
#%%
plt.figure('Eigenvalue errors_2')
i_g = 21
plt.ylabel('|E_a - E_num|/E_a in %')
plt.xlabel('Mode')
plt.plot(np.linspace(0,len(E)-1,len(E), dtype= int)[0:i_g], np.abs(E_analytic[0:i_g]-E[0:i_g])/E_analytic[0:i_g]*100, 'rx', markersize = '8.5')
plt.xticks(np.linspace(0,len(E)-1,int(len(E)*2), dtype= int)[0:i_g*2])
#plt.savefig('assgn4_abse_percent.pdf', dpi=1000,bbox_inches='tight')
#%%
