# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:24:19 2020

@author: CL
"""
import chem_RW as crw
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp

def multisim(func, L, N, t, steps, sims):
    data0 = func(L,N,t,steps)
    for j in range (0, sims):
        data = func(L, N, t,steps)
        data0 = np.vstack([data0, data])
    
    avg = np.mean(data0, axis = 0)
    
    return data0, avg


def datainv(array):
    data = []
    for i in array:
        if i > 0:
            data += [1/i]
    return data


oneDrun = multisim(crw.one_D,729, 650, 1000, 1, 10)[1]
t1 = (np.arange(0,len(oneDrun)))  
plt.plot(t1, oneDrun)
plt.title("1D Lattice: [A] vs Time") 
plt.figure()

twoDrun = multisim(crw.two_D,27, 650, 250, 1, 10)[1]
t2 = (np.arange(0,len(twoDrun)))
plt.plot(t2, twoDrun)
plt.title("2D Lattice: [A] vs Time")
plt.figure()

threeDrun = multisim(crw.three_D, 9, 650, 100, 1, 10)[1]
t3 = (np.arange(0,len(threeDrun)))
plt.plot(t3, threeDrun)
plt.title("3D Lattice: [A] vs Time")
plt.figure()

x1, y1 = t1[0:700], datainv(oneDrun)[0:700]
x2, y2 = t2[0:150], datainv(twoDrun)[0:150]
x3, y3 = t3[0:70], datainv(threeDrun)[0:70]

m1, b1, r1 = sp.linregress(x1, y1)[0:3]
m2, b2, r2 = sp.linregress(x2, y2)[0:3]
m3, b3, r3 = sp.linregress(x3, y3) [0:3]
  
print("Fit 1")
print (r1**2)
plt.title("1D Lattice: 1/[A] vs Time")
oneDinv = datainv(oneDrun)
t1a = (np.arange(0,len(oneDinv))) 
plt.plot(t1a[0:700], oneDinv[0:700], 'bo', markersize = 2)
plt.plot(t1a[0:700], m1*t1a[0:700] + b1, 'r')
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()

print("Fit 2")
print (r2**2)
plt.title("2D Lattice: 1/[A] vs Time")
twoDinv = datainv(twoDrun)
t2a = (np.arange(0,len(twoDinv)))
plt.plot(t2a[0:150], twoDinv[0:150], 'bo', markersize = 2) 
plt.plot(t2a[0:150], m2*t2a[0:150] + b2, 'r')
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()

print("Fit 3")
print (r3**2)
plt.title("3D Lattice: 1/[A] vs Time")
threeDinv = datainv(threeDrun)
t3a = (np.arange(0,len(threeDinv))) 
plt.plot(t3a[0:70], threeDinv[0:70], 'bo', markersize = 3)
plt.plot(t3a[0:70], m3*t3a[0:70] + b3, 'r')
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()
