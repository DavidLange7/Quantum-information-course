#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 18:42:00 2022

@author: david
"""
import numpy as np

class assgn6:
    
    def __init__(self, N, d, sep_check, psi_coeff = None, debug = True):
        self.psi_coeff = psi_coeff
        self.d = d
        self.N = N            
        self.sep_check = sep_check
        self.debug = debug
        
    def init_wf(self):
        if self.sep_check == True:
            shape = self.d*self.N
        else:
            shape = self.d**(self.N)
        
        if self.psi_coeff is None:
            psi = np.array(np.random.uniform(-1, 1, shape) + 1j * np.random.uniform(-1, 1, shape))
        else:
            psi = self.psi_coeff
        if self.debug == True:
            print('---print raw psi---')
            print(psi)

        nm = np.dot(psi,psi.conj())
        nm = np.real(nm)
        psi = psi/np.sqrt(nm)
        
        if self.debug == True:
            print('--------------print psi normalized--------------')
            print(psi)
            print('normalized?')
            norm_i = np.sum(np.dot(psi,psi.conj()))
            
            print(norm_i)
        return(psi)
    
    def density_mat(self):
        psi_t = self.init_wf()
        dens = np.array([[psi_t[j]*np.conj(psi_t[i]) for i in range(len(psi_t))] for j in range(len(psi_t))])
        if self.debug == True:
            print('---------------------------')
            print(dens)
            print('-----check if p^2-p=0------')
            print(np.dot(dens,dens)-dens)
            print('-----check if tr(p^2) = tr(p) = 1------')
            print(np.trace(np.dot(dens,dens)), np.trace(dens))
        return(dens)
            
            
    def density_mat_red(self, sys):
        self.debug = False
        dens = self.density_mat()
        r_mat = np.zeros([self.d**(self.N-1), self.d**(self.N-1)], dtype=complex)
        
        for i in range(self.d**(self.N-1)): #rows
            for j in range(self.d**(self.N-1)): #columns

                tmp = 0 + 1j*0
                for n in range(self.d):

                    if sys == 'L':
                        alph_1 = 2*i + n
                        alph_2 = 2*j + n
                        tmp += dens[alph_1, alph_2]
                        if self.debug == 'True':
                            print([alph_1, alph_2])
                        
                    if sys == 'R':
                        alph_1 = 2*n + i
                        alph_2 = 2*n + j
                        tmp += dens[alph_1, alph_2]
                        if self.debug == 'True':
                            print([alph_1, alph_2])
                r_mat[i, j] = tmp
        print(r_mat)