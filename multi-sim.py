# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 12:25:50 2020

@author: CL
"""
import chem_RW as crw
import numpy as np
import matplotlib.pyplot as plt

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

######## 1D Simulations
data1 = crw.one_D(4096, 250, 12000, 500)

for i in range (0,10):
    data = crw.one_D(4096, 250, 12000, 500)
    data1 = np.vstack([data1, data])

print(data1)

avg1 = np.mean(data1, axis = 0)


t1 = 500*(np.arange(0,len(avg1)))  
sigma1 = stdev(data1, avg1, 10)
plt.errorbar(t1, avg1, sigma1)
plt.figure()
######## 2D Simulations

data2 = crw.two_D(64, 250, 2500, 100)

for k in range (0,10):
    data = crw.two_D(64, 250, 2500, 100)
    data2 = np.vstack([data2, data])

print(data2)
avg2 = np.mean(data2, axis = 0)

t2 = 100*(np.arange(0,len(avg2)))
sigma2 = stdev(data2, avg2, 10)
plt.errorbar(t2, avg2, sigma2)
plt.figure()

######## 3D Simulations

data3 = crw.three_D(16, 250, 1500, 50)

for j in range (0,10):
    data = crw.three_D(16, 250, 1500, 50)
    data3 = np.vstack([data3, data])

print(data3)
avg3 = np.mean(data3, axis = 0)

t3 = 50*(np.arange(0,len(avg3)))
sigma3 = stdev(data3, avg3, 10)
plt.errorbar(t3, avg3, sigma3)
plt.figure()