#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 15:10:39 2022

Purpose of the code: Implements Assignment 5 of the quantum information course autumn 2022, split operator method to solve time-dependent SE.

@author: david
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from matplotlib import cm
import datetime
from scipy import optimize

plt.style.use('ggplot')
#%%

class tdse:
    def __init__(self, T_max = 430, xmin = -10, xmax = 10, n_mesh = 1000, n_t = 1000, imag_time = False):
        '''
        Parameters
        ----------
        T_max : TYPE, float
            is the maximum time of propagation
        dt : TYPE, float
            time step
        xmin, xmax : TYPE, float
            edges of the space considered
        n_mesh ; TYPE, integer
            number of points of discretized space
        n_t : TYPE, integer
            number of points of discretized time
        '''
        self.T_max = T_max
        self.dt = T_max/n_t
        self.n_t = n_t
        self.xmin = xmin
        self.xmax = xmax
        self.n_mesh = n_mesh
        
        self.h = (xmax-xmin)/n_mesh
        self.xsp = np.array([xmin + self.h*i for i in range(0,n_mesh)])
        self.ksp = np.fft.fftfreq(n_mesh, self.h)*2*np.pi
        
        if imag_time == True:
            self.var_temp = 1
        else:
            self.var_temp = 1j
        pass

        
    def __inputf(self, x, a = 1):
        '''
        Input wave packet, for this code the ground state of the harmonic oscillator
        '''
        return (np.array([(a/np.pi)**(1/4)*np.exp(-np.square(x[i])/2) for i in range(0,len(x))]))
    
    def __U_k(self, psi_k,k,dt):
        '''
        The kinetic term evolution of the hamiltonian, in fourier space
        '''
        return(psi_k*np.exp(-self.var_temp*k**2/2*dt))
    
    def __U_x(self, psi, x, dt, T, t):
        '''
        The potential term evolution of the hamiltonian, in real space
        '''
        return(psi*np.exp(-self.var_temp*(x-t/T)**2/2*dt/2))

    def comp(self, plot1 = False, plot_video = False, print_norms = False):
        '''
        The split operator computation, plus some plotting/printing ptions
        print_norms = True can be nice to check if input/output wavefunctions are normalized correctly as a checkpoint..
        '''
        psi = np.zeros((self.n_t,len(self.xsp)), dtype = complex)
        
        psi[0,:] = self.__inputf(self.xsp)
        
        nm = np.dot(psi[0,:],psi[0,:].conj())*self.h
        nm = np.real(nm)
        psi[0,:] = psi[0,:]/np.sqrt(nm)

        
        if plot1 == True:
            fig2, (ax21, ax22) = plt.subplots(2)
            ax21.plot(self.xsp, np.abs(psi[0,:])**2)
            ax21.plot(self.xsp, (self.xsp)**2/2)
            ax21.set_ylim([-1, 1])
        
        norm_i = np.sum(np.dot(psi[0,:],psi[0,:].conj())*self.h)
        if print_norms == True:
            print('Is the initial wavepacket normalized:', norm_i)
        
        t = 0
        if plot_video == True:
            fig, (ax1) = plt.subplots(1)
            res = np.abs(psi[0,:])**2
            line1, = ax1.plot(self.xsp,res,'black')
            line2, = ax1.plot(self.xsp,(self.xsp-self.dt/self.T_max)**2/2,'g-')        
            plt.xlabel('Time')
            plt.ylabel('|\u03C8|²')

        for i in range(1, self.n_t):
            
            nm = np.dot(psi[i-1,:],psi[i-1,:].conj())*self.h
            nm = np.real(nm)
            psi[i-1,:] = psi[i-1,:]/np.sqrt(nm)
            
            psi[i,:] = self.__U_x(psi[i-1,:], self.xsp, self.dt, self.T_max, self.dt*i)
            
            psi_k = (np.fft.fft(psi[i,:]))
            
            psi_ki = self.__U_k(psi_k, self.ksp, self.dt)
            psi[i,:] = np.fft.ifft((psi_ki))
            
            nm = np.dot(psi[i,:],psi[i,:].conj())*self.h
            nm = np.real(nm)
            psi[i,:] = psi[i,:]/np.sqrt(nm)
            
            psi[i,:] = self.__U_x(psi[i,:],self.xsp,self.dt,self.T_max, self.dt*i)
            
            
            if plot_video == True:
                res = np.abs(psi[i,:])**2
                line1.set_ydata(res)
                line2.set_ydata((self.xsp-t/self.T_max)**2/2)
                ax1.set_ylim([0, 1])
                ax1.set_xlim([-5, 5])
                ph = plt.fill_between(self.xsp, (self.xsp-t/self.T_max)**2/2, color='g',alpha=.5)
                fig.canvas.draw()
                fig.canvas.flush_events()
                plt.pause(0.0005)
                '''
                #This part is only used to create multiple pngs that can be put together in a gif using imageio ..
                
                plt.savefig('%d.png'% (i))
                file = open("/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment5/pics.txt", "a")
                file.write('"')
                file.write('%d.png'% (i))
                file.write('"')
                file.write(',')
                file.close()
                '''
                if i != math.ceil(self.T_max/self.dt)-1:
                    ph.remove()
            t += self.dt
            
        norm_f = np.sum(np.dot(psi[i,:],psi[i,:].conj())*self.h)
        
        if print_norms == True:
            print('Is the final wavepacket normalized:', norm_f)
        res = np.abs(psi[-1,:])**2
        
        if plot1 == True:
            ax22.plot(self.xsp, res, color = 'gray')
            ax22.set_ylim([-0.1, np.max(res)])
        return(psi)
    
    def contour(self):
        '''
        Plots the whole time evolution of the wavefunction in a nice contour plot
        '''
        psi = self.comp(plot1 = False, plot_video = False)
        
        fig3, ax3 = plt.subplots()
        cs = ax3.contourf(self.xsp, np.linspace(0,self.T_max, self.n_t), np.abs(psi**2), levels=60, cmap=cm.gnuplot, extend='max')
        cbar = fig3.colorbar(cs)#, ticks=[if you want more customization])
        plt.axvline(x=1, linestyle = 'dashed', color = 'gray', alpha = 0.6)
        ax3.set_xlim(-5,5)
        ax3.set_xlabel('x', fontsize = 11)
        ax3.set_ylabel('time', fontsize = 11)

        cbar.ax.set_ylabel('|\u03C8|²', fontsize = 11)
        #plt.savefig('contour_assgn5_qi.pdf', dpi=1500,bbox_inches='tight')
        pass
        
    def deviation(self, T_man = 0, plot = False, frequency_stuff = False):
        '''
        Basically takes one T_max value as input and computes <x(t)>, can plot if wanted..
        returns the max deviation between <x(t)> and expected behaviour which is x = t/T
        
        Also, frequency_stuff = True gives an additional option to check the frequency change of the oscillations with increasing T_max ...
        '''
        if T_man != 0:
            self.T_max = T_man
            self.dt = self.T_max/self.n_t
            
        psi = self.comp(plot1 = False, plot_video = False)
        psi_sq = np.abs(psi**2)
        
        expc_x = np.array([np.dot(psi_sq[i,:], self.xsp*self.h) for i in range(self.n_t)])
        
        
        if frequency_stuff == True:
            def g(x, A, w):
                return(A*np.sin(w*x))
            popt, pcov = optimize.curve_fit(g, np.linspace(0,self.T_max, self.n_t), expc_x-np.linspace(0,1, self.n_t))
            plt.plot(np.linspace(0,self.T_max, self.n_t), expc_x - np.linspace(0,1, self.n_t), label = 'data')
            plt.plot(np.linspace(0,self.T_max, self.n_t), g(np.linspace(0,self.T_max, self.n_t), popt[0], popt[1]), label = 'fit')
            plt.legend()
            #remember to swap return for return float(popt[1]) to get frequency returned for meanstuff() function
        
        if plot == True:
            plt.plot(np.linspace(0,self.T_max, self.n_t), np.linspace(0,1, self.n_t), label = 'expected')
            plt.plot(np.linspace(0,self.T_max, self.n_t), expc_x, label = '<x(t)>')
            plt.xlabel('time')
            plt.ylabel('space')
            plt.legend()
            #plt.savefig('deviat_assgn5_qi_time250.pdf', dpi=1500,bbox_inches='tight')
            
        max_dev = np.max(np.abs(expc_x-np.linspace(0,1, self.n_t)))
        return max_dev
    
    
    def meanstuff(self):
        '''
        Takes self.deviation() and runs it for a space of T_max, plus plot and fit on the result
        '''
        tsp = np.linspace(10,400,50)
        max_devs = (np.array([self.deviation(T_man = tsp[i], plot = False) for i in range(len(tsp))]))
        
        plt.plot(tsp, max_devs, 'rx')
        def f(x, a, b):
            return(a*x**b)
        popt, pcov = optimize.curve_fit(f, tsp, max_devs)
        plt.plot(tsp, f(tsp, popt[0], popt[1]), color = 'gray')
        plt.xlabel('T_max')
        plt.ylabel('Max deviation')
        
        print(popt[0], popt[1])
        
        #plt.savefig('meanstuff_assgn5_qi.pdf', dpi=1500,bbox_inches='tight')
        pass

#%%
#this is for making a gif out of some png files ..
'''
import imageio
import os

os.chdir('/home/david/Courses_padova/Last semester FUN stuff/QuantumComp/Assignment5')
with imageio.get_writer('video_time.gif', mode='I') as writer:
    for filename in [...]:
        image = imageio.imread(filename)
        writer.append_data(image)
'''