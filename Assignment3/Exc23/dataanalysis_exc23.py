#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 13:27:36 2022

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
import pylab as plb

plt.style.use('ggplot')
plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (16,10)
#%%

spc = np.loadtxt("/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment3/normspacing.txt")
spc = spc.flatten()
spc = [x for x in spc if x <= 4]
#%%
#S = [spc[i+1]-spc[i] for i in range(0,len(spc)-1)]
#avg = sum(S)/6999
#S = S/avg
#%%
plt.figure(2)
hist = plt.hist((spc), bins = 500, density = 'True', stacked = 'True')

#problem could actually be that you have to increase sensitivity, by increasing number of bins ... 
def f(s, a, b, al, be):
    return (a*s**al*np.exp(b*s**be))

res = hist[0]/sum(hist[0])
#%%
i_temp = 300

#plt.plot(hist[1][:i_temp], res[:i_temp])
popt, pcov = cf(f, hist[1][:i_temp], res[:i_temp])

plt.figure(1)
plt.plot(hist[1][:i_temp], res[:i_temp], 'bx', alpha = 0.6)
plt.ylabel('P(s)')
plt.xlabel('s')
plt.plot(hist[1][:i_temp], f(hist[1][:i_temp], popt[0], popt[1], popt[2], popt[3]), linewidth = 2)
#plt.savefig('assgn3_exc3_herm.pdf', dpi=1000,bbox_inches='tight')
plt.show()