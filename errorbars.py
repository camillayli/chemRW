# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 12:25:50 2020

@author: CL
"""
import chem_RW as crw
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp


def stdev (indiv, avgs, sims):
    sigmas = []
    print(avgs.size)
    print(sims)
    for m in range (0, avgs.size): 
        sumsi = []
        
        for n in range (0, sims):
            sumsi += [(indiv[n][m]-avgs[m])**2]
            
        sumstotal = np.sum(sumsi)
        sdev = sumstotal/(sims - 1)
        
        sigmas += [sdev]
        
    return sigmas

 
def multisim(func, L, N, t, steps, sims):
    data0 = func(L,N,t,steps)
    for j in range (0, sims):
        data = func(L, N, t,steps)
        data0 = np.vstack([data0, data])
    
    avg = np.mean(data0, axis = 0)
    
    return data0, avg




data1, oneDavg = multisim(crw.one_D,10, 8, 10, 1, 10)


sigma1 = stdev(data1, oneDavg, 10)
t1a = np.arange(0, len(oneDavg))


plt.errorbar(t1a, oneDavg, sigma1)
plt.figure()
