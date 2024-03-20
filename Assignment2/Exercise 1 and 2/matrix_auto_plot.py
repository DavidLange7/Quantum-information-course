#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')
import pylab as plb
plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (16,10)


A_1 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assigment1/one.txt')
A_2 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assigment1/two.txt')
A_3 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assigment1/three.txt')

n_0 = [1]
n_1 = np.linspace(15,50,36)
n_2 = 50*np.linspace(2, 30, 29)
n = np.concatenate((n_0, n_1, n_2))
#%%
plt.yscale('log')
plt.xscale('log')

#plt.title('-O0 - optimization flag')
plt.plot(n, A_1, label = 'rowwise', color = 'red', linewidth = '2')
plt.plot(n, A_2, label = 'columnwise', color = 'blue', linewidth = '2')
plt.plot(n, A_3, label = 'matmul', color = 'green', linewidth = '2')
plt.grid()
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.legend()

'''
#axes = plt.axes([left, bottom, width, height])
axes = plt.axes([0.528, 0.2, 0.35, 0.25])
axes.plot(n, A_1 , color = 'red')
#axes.plot(n, A_1, 'rx')
axes.plot(n, A_2 , color = 'blue')
axes.plot(n, A_3 , color = 'green')
plt.yscale('log')
#plt.yticks([0,5,10,15])
plt.xlim([15, 50])
plt.ylim([10**(-6), 10**(-4)])
plt.grid()
'''
#plt.savefig('qi_ass1_00_funn_optimization.pdf', dpi=1000,bbox_inches='tight')
#%%
def f(x,m,b):
    return (m*x + b)

#first value is kind of weird
n_temp = n[1:]
A_1_temp = A_1[1:]
A_3_temp = A_3[1:]

res1 = np.polyfit(np.log(n_temp), np.log(A_1_temp), 1)
res2 = np.polyfit(np.log(n_temp), np.log(A_3_temp), 1)

plt.figure('fit1')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.plot(n_temp, A_1_temp, 'rx', label = 'rowwise data')
plt.plot(n_temp, np.exp(f(np.log(n_temp), res1[0], res1[1])))
plt.grid()
plt.legend()
#plt.savefig('sfdghjk.pdf', dpi=1000,bbox_inches='tight')

plt.figure('fit2')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.plot(n_temp, A_3_temp, 'gx', label = 'matmul data')
plt.plot(n_temp, np.exp(f(np.log(n_temp), res2[0], res2[1])))
plt.grid()
plt.legend()
#plt.savefig('sfdghjk.pdf', dpi=1000,bbox_inches='tight')