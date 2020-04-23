# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:09:47 2020

@author: CL
"""

import numpy as np
import random as rd
import numpy.random as nprnd
import matplotlib.pyplot as plt



def remove_dup(array):
    """
    Function to remove elements that are duplicates in a numpy array. 
    Returns a new numpy array with those elements removed and the number of 
    remaining molecules.
    If there are no duplicate elements in the input array, returns the input array.
    If there are no more elements in the array, returns empty array.
    """
    if (array.size > 0):
        coparray, dupe = np.unique(array, axis = 0, return_counts = True)
        print(dupe)
        dupe2 = []
            
            
        for i in range(0, dupe.size):
            if dupe[i] > 1:
                dupe2 += [i]
                    
        coparray = np.delete(coparray, dupe2, 0)
        if (len(coparray.shape) == 1):
            size = coparray.size
        else:    
            size = coparray[:,0].size
           
        
    else:
        size = 0
        print(0)
        a = []
        coparray = np.array(a)
        
    return coparray, size


def one_D(L, N, trials, steps):
    """
    1-D lattice  
    Assumptions: Molecule must move by one cell in each trial;
        they won't stop moving just because they are at the boundary.
        Molecules annihilate only when they occupy the same space.
    Boundary conditions for 1-D lattice: 
        if index is 0 or L-1, they move away from boundary.
    
    """
    coord = [] #list to store coordinates of molecules
    data = [N/L] #concentration
    
    for i in range (0, N):
        x = rd.randrange(0, L) #generate random integers - starting positions - x coordinate
        while x in coord: #make sure there are no overlapping initial positions
            x = rd.randrange(0, L)

        coord += [x]
    
    coordcopy = np.array(coord)

    for n in range (0, int(trials/steps)): 
        for j in range (0, steps):
            
            
            for k in range (0, coordcopy.size):
                chance = rd.uniform(0,1)  #random float from 0 to 1
            
                #implements boundary conditions
                if ((chance < 0.5) and (coordcopy[k]!= 0)) or (coordcopy[k] == L - 1): 
                    coordcopy[k] -= 1
                    
                else:
                    coordcopy[k] += 1
                
            coordcopy, count = remove_dup(coordcopy)
         
        data += [count/L]
        
    return data  

 
def two_D(L, N, trials, steps):
    coord = []
    data = [N/(L**2)]
    
    for i in range (0, N):
        x = rd.randrange(0, L)  #generate random integers - starting positions
        y = rd.randrange(0, L)  #x and y coordinates
        while [x,y] in coord:
            x = rd.randrange(0, L)  
            y = rd.randrange(0, L)
            
        coord += [[x, y]]
    
    coordcopy = np.array(coord)
    
    for n in range (0, int(trials/steps)):
        for j in range (0, steps):
            
            if (coordcopy.size != 0):
                
                for k in range (0, coordcopy[:, 0].size):
                    
                    chance = (nprnd.random(3))  # 3 random floats from 0 to 1
                    #print(chance)
                    #scenario where molecule moves in both x and y directions
                    if (chance[2] < 0.5):
                        indices = [0, 1]
                    #scenario where molecule moves in only the x direction
                    elif (chance[2] < 0.75): 
                        indices = [0]
                    #scenario where molecule moves in only the y direction   
                    else:
                        indices = [1]
                        
                    for m in indices: 
                        if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                            coordcopy[:, m][k] -= 1
                                
                        else: 
                            coordcopy[:, m][k] += 1
                         
            #implements the annihilation by updating coordcopy
            #removes coordinates that are duplicates 
            coordcopy, count = remove_dup(coordcopy)
      
        data += [count/(L**2)]
        
    return data
        

def three_D(L, N, trials, steps):
    coord = []
    data = [N/(L**3)]
    
    for i in range (0, N):
        x = rd.randrange(0, L)  #generate random starting positions
        y = rd.randrange(0, L)  #x, y, and z coordinates
        z = rd.randrange(0, L)
        while [x,y,z] in coord:
            x = rd.randrange(0, L)  
            y = rd.randrange(0, L)  
            z = rd.randrange(0, L)
   
        coord += [[x, y, z]]
    
    coordcopy = np.array(coord)
    
    
    for n in range (0, int(trials/steps)):
        for j in range (0, steps):
            
            if (coordcopy.size != 0):
                for k in range (0, coordcopy[:, 0].size):
                    
                    chance = (nprnd.random(4))  # 4 random floats from 0 to 1
                    #print(chance)
                    #scenario where molecule moves in x, y, and z directions
                    if (chance[3] < (4.0/13)):
                        indices = [0,1,2]
                    #scenario where molecule moves in only the x and y directions
                    elif (chance[3] < (6.0/13)): 
                        indices = [0,1]
                    #scenario where molecule moves in only the x and z directions   
                    elif (chance[3] < (8.0/13)):
                        indices = [0,2]
                    #scenario where molecule moves in only the y and z directions
                    elif (chance[3] < (10.0/13)):
                        indices = [1,2]
                    #scenario where molecule moves in only the x direction
                    elif (chance[3] < (11.0/13)):
                        indices = [0]
                    #scenario where molecule moves in only the y direction
                    elif (chance[3] < (12.0/13)):
                        indices = [1]
                    #scenario where molecule moves in only the z direction
                    else:
                        indices = [2]
                    
                    for m in indices: 
                        if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                            coordcopy[:, m][k] -= 1
                                
                        else: 
                            coordcopy[:, m][k] += 1
            
            #implements the annihilation by updating coordcopy
            #removes coordinates that are duplicates 
            coordcopy, count = remove_dup(coordcopy)
         
        data += [count/(L**3)]
        
    return data  

"""                        
oneDrun = one_D(4096, 250, 12000, 500)
t = 500*(np.arange(0,len(oneDrun)))  
plt.plot(t, oneDrun)
plt.title("1D Lattice: 1/[A] vs Time") 
plt.figure()

twoDrun = two_D(64, 250, 2500, 100)
t = 100*(np.arange(0,len(twoDrun)))
plt.plot(t, twoDrun)
plt.title("2D Lattice: 1/[A] vs Time")
plt.figure()

threeDrun = three_D(16, 250, 1500, 50)
t = 50*(np.arange(0,len(threeDrun)))
plt.plot(t, threeDrun)
plt.title("3D Lattice: 1/[A] vs Time")
plt.figure()
"""