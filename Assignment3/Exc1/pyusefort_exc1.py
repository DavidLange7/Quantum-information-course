#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 08:39:02 2022

@author: david
"""
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
plt.style.use('default')
import pylab as plb
plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (16,10)

#tmp1 = np.linspace(1,100,100) #RUN THIS AGAIN IF YOU HAVE TIME AND PUT INTO THE OVERLEAF
tmp2 = np.linspace(250,3000,70)
#n = np.concatenate((tmp1, tmp2))
n = tmp2.astype(int)
#%%
call(["bash", "reset.sh"])
for i in n:
    with open('settings.txt', 'w') as f:
        f.write('%.2d'%(i))
    call(["gfortran","-o","out","matrix.f90"])
    call(["/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment3/out"])
    #%%
A_1 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment3/one.txt')
A_2 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment3/two.txt')
A_3 = np.loadtxt('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment3/three.txt')

#plt.yscale('log')
#plt.xscale('log')
def f(x,m,b):
    return (m*x + b)

def g(x,m,b):
    return(m*x**b)

#first value is kind of weird
n_temp = n
A_1_temp = A_1
A_3_temp = A_3

popt1, pcov1 = cf(f, np.log(n), np.log(A_1))
popt2, pcov2 = cf(f, np.log(n), np.log(A_3))
popt3, pcov3 = cf(f, np.log(n), np.log(A_2))


plt.figure('fit1')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.plot(n_temp, A_1_temp, 'gx', label = 'rowwise data')
plt.plot(n_temp, np.exp(f(np.log(n_temp), popt1[0], popt1[1])), color = 'green')
plt.grid()
plt.legend()
#plt.savefig('sfdghjk.pdf', dpi=1000,bbox_inches='tight')

#plt.figure('fit2')
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.plot(n_temp, A_3_temp, 'rx', label = 'matmul data')
plt.plot(n_temp, np.exp(f(np.log(n_temp), popt2[0], popt2[1])), color = 'red')
plt.grid()
plt.legend()
#plt.savefig('sfdghjk.pdf', dpi=1000,bbox_inches='tight')

#plt.figure('fit3')
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Matrix dimension')
plt.ylabel('Computation time in [s]')
plt.plot(n_temp, A_2, 'bx', label = 'columnwise data')
plt.plot(n_temp, np.exp(f(np.log(n_temp), popt3[0], popt3[1])), color = 'blue')
plt.grid()
plt.legend()
#plt.savefig('assgn3_exc1.pdf', dpi=1000,bbox_inches='tight')