# -*- coding: utf-8 -*-
"""
Created on Fri Apr  10 13:53:59 2020

@author: Camilla Li
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
    while True:
        try:
            coparray, dupe = np.unique(array, axis = 0, return_counts = True)
            print(dupe)
            dupe2 = []
               
            for i in range(0, dupe.size):
                if dupe[i] > 1:
                    dupe2 += [i]
                    
            coparray = np.delete(coparray, dupe2, 0)
            size = coparray[:,0].size
            break
        
        except ValueError:
            size = 0
            coparray = np.array()
        
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
    array1 = np.zeros(L) #lattice where simulation occurs
    for i in range (0, N):
        x = rd.randrange(0, L) #generate random integers - starting positions - x coordinate
        while (array1[x] != 0): #make sure there are no overlapping initial positions
            x = rd.randrange(0, L)
        
        array1[x] = 1
        coord += [x]
    
    print(array1)
    print(coord)
    
    data = [N]
    count = N
    #for n in range (0, steps): 
    for j in range (0, trials):
        #temparray = np.copy(array1) 
        
        
        for k in coord:
            chance = rd.uniform(0,1)  #random float from 0 to 1
        
            #implements boundary conditions
            if ((chance < 0.5) and (k!= 0)) or (k == L - 1): 
                array1[k - 1] += 1 
                coord[coord.index(k)] -= 1
                
            else:
                array1[k + 1] += 1
                coord[coord.index(k)] += 1
            
            #to remove from its previous position
            array1[k] -= 1
        
        #implements the annihilation - only when molecules occupy the same space
        #this occurs only after the for loop above is complete, and all the molecules have moved to their new positions
        
        for m in range (0, L):
            if array1[m] == 2:
                array1[m] = 0 
                coord.remove(m) #deletes coordinate of annihilated molecule
                coord.remove(m) #each interation only removes the first instance
                count -= 2 #updates counter
         
        print(array1)        
        print(coord)
        
        data += [count]
        
    return data

 
def two_D(L, N, trials, steps):
    array2 = np.zeros((L,L))
    coord = []
    data = [N]
    
    for i in range (0, N):
        x = rd.randrange(0, L)  #generate random integers - starting positions
        y = rd.randrange(0, L)  #x and y coordinates
        while [x,y] in coord:
            x = rd.randrange(0, L)  
            y = rd.randrange(0, L)
            
        
        array2[x][y] = 1
        coord += [[x, y]]
    
    coordcopy = np.array(coord)
    print(array2)
    
    print(coordcopy)
    
    #for n in range (0, steps):
    for j in range (0, trials):
        #temparray = np.copy(array2)
        
        for k in range (0, coordcopy[:, 0].size):
            
            chance = (nprnd.random(3))  # 3 random floats from 0 to 1
            print(chance)
            #scenario where molecule moves in both x and y directions
            if (chance[2] < 0.5):
                index1 = 0
                index2 = 2
            #scenario where molecule moves in only the x direction
            elif (chance[2] < 0.75): 
                index1 = 0
                index2 = 1
            #scenario where molecule moves in only the y direction   
            else:
                index1 = 1
                index2 = 2
                
            for m in range (index1, index2): 
                if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                    coordcopy[:, m][k] -= 1
                        
                else: 
                    coordcopy[:, m][k] += 1
        
        print(coordcopy)                
        #implements the annihilation by updating coordcopy
        #removes coordinates that are duplicates 
        coordcopy, count = remove_dup(coordcopy)
        
        array2 = np.zeros((L,L))
        for p in range (0, coordcopy[:,0].size):
            array2[coordcopy[p][0]][coordcopy[p][1]] = 1
          
        print(array2)
        print(coordcopy)  
        data += [count]
        
    return data
        

def three_D(L, N, trials, steps):
    coord = []
    data = [N]
    
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
    
    print(coordcopy)
    
    #for n in range (0, steps):
    for j in range (0, trials):
        #temparray = np.copy(array2)
        
        for k in range (0, coordcopy[:, 0].size):
            
            chance = (nprnd.random(4))  # 4 random floats from 0 to 1
            print(chance)
            #scenario where molecule moves in x, y and z directions
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
        
        print(coordcopy)  
        data += [count]
        
    return data  
                        
oneDrun = one_D(500, 250, 100, 10)
t = np.arange(0,len(oneDrun))  
plt.plot(t, oneDrun)
plt.title("1D Lattice: Concentration vs Time") 
plt.figure()

twoDrun = two_D(80, 250, 500, 10)
t = np.arange(0,len(twoDrun))
plt.plot(t, twoDrun)
plt.title("2D Lattice: Concentration vs Time")
plt.figure()

threeDrun = three_D(20, 250, 800, 10)
t = np.arange(0,len(threeDrun))
plt.plot(t, threeDrun)
plt.title("3D Lattice: Concentration vs Time")
plt.figure()
