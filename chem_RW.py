# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:09:47 2020

@authors: Camilla Li, Noah Walton, Matthew Willard
"""
#Importing of built-in libraries
import numpy as np
import random as rd
import numpy.random as nprnd
import matplotlib.pyplot as plt
import scipy.stats as sp



######## Part 1. Defining Functions for Simulating Chemical Annihilation in 1, 2, and 3 Dimensions
def one_D(L, N, t, steps):
    """
    1-D simulation  
    Assumptions: Particles must move by one cell in each trial;
        they won't stop moving just because they are at the boundary.
        When pairs of particles occupy the same space, they annihilate.
    Boundary conditions for 1-D lattice: 
        If index is 0 or L-1, when the particle would have "crossed over" the 
        boundary, they instead move away from the boundary.
    Inputs: 
        L = size of grid
        N = number of initial walkers
        t = time steps for which we will run the simulation
        Steps: number of time steps between which we will make measurements of 
        concentration
    
    Outputs measurements of concentration (of walkers) over time. 
    
    """
    coord = [] #list to store coordinates of partocles 
    data = [N/L] #initial concentration of random walkers
    
    for i in range (0, N):
        x = rd.randrange(0, L) #generate random integers - starting positions - x coordinate
        while x in coord: #make sure there are no overlapping initial positions
            x = rd.randrange(0, L)

        coord += [x] #concatenates unique coordinate into list
    
    coordcopy = np.array(coord) #converts list to a numpy array, so array operations can be performed

    for n in range (0, int(t/steps)): 
        for j in range (0, steps):
            
            for k in range (0, coordcopy.size): 
                #loops through the coordinates of the random walkers
                #if there are no more walkers left, the loop is not executed
                
                chance = rd.uniform(0,1)  #random float from 0 to 1
            
                #implements boundary conditions
                #coordinate can increment or decrement depending on random number 
                #or its location if it is bordering the boundary
                if ((chance < 0.5) and (coordcopy[k]!= 0)) or (coordcopy[k] == L - 1): 
                    coordcopy[k] -= 1
                    
                else:
                    coordcopy[k] += 1
                    
            #implements the annihilation by updating coordcopy
            #removes pairs of coordinates that are duplicates - ie. pair of random walkers in the same space     
            #the function remove_dup is defined later in this file
            coordcopy, count = remove_dup(coordcopy)
         
        data += [count/L] #records the concentration of walkers after the number of 
                            #time steps specified by the variable "steps". 
        
    return data

 
def two_D(L, N, t, steps):
    """
    2-D simulation  
    Assumptions: Particles must move by one cell in each trial;
        they won't stop moving just because they are at the boundary.
        When pairs of particles occupy the same space, they annihilate.
    Boundary conditions for 2-D lattice: 
        If index for any coordinate is 0 or L-1, when the particle would have 
        "crossed over" the boundary, they instead move away from the boundary.
    Inputs: 
        L = length of LxL grid
        N = number of initial walkers
        t = time steps for which we will run the simulation
        Steps: number of time steps between which we will make measurements of 
        concentration
    
    Outputs measurements of concentration (of walkers) over time. 
    """
    coord = []
    data = [N/(L**2)] #initial concentration of walkers
    
    for i in range (0, N):
        x = rd.randrange(0, L)  #generate random integers - starting positions
        y = rd.randrange(0, L)  #x and y coordinates
        while [x,y] in coord:  
            x = rd.randrange(0, L)  
            y = rd.randrange(0, L)
        #if there is already a random walker occupying a certain location, 
        #new coordinates for x and y are generated until the walker is given
        #a unique starting position
            
        coord += [[x, y]] #concatenates unique coordinate into list
    
    coordcopy = np.array(coord) 
    #list is converted to np array so array operations can be performed
    
    for n in range (0, int(t/steps)):
        for j in range (0, steps):
            
            if (coordcopy.size != 0):
                
                for k in range (0, coordcopy[:, 0].size): 
                #loops through the coordinates of the random walkers
                #if there are no more walkers left, the loop is not executed
                
                    chance = (nprnd.random(3))  
                    #generates 3 random floating point numbers from 0 to 1
                    
                    #scenario where particle moves in both x and y directions
                    if (chance[2] < 0.5):
                        indices = [0, 1]
                    #scenario where particle moves in only the x direction
                    elif (chance[2] < 0.75): 
                        indices = [0]
                    #scenario where particle moves in only the y direction   
                    else:
                        indices = [1]
                        
                    #implements boundary conditions
                    #loops through the values in the list "indices," so that for 
                    #the coordinate(s) along which the particle can move (determined
                    #by a random number), each coordinate(s) can increment or decrement 
                    #depending on separate random number or its location if it is bordering the boundary
                    #Each coordinate is incremented/decremented independently from the other.
                    #This method allows for scalability when changing the dimensions
                    for m in indices: 
                        if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                            coordcopy[:, m][k] -= 1
                                
                        else: 
                            coordcopy[:, m][k] += 1
                         
            #implements the annihilation by updating coordcopy
            #removes pairs of coordinates that are duplicates - ie. pairs of random walkers in the same space 
            coordcopy, count = remove_dup(coordcopy)
      
        data += [count/(L**2)] #records the concentration of walkers after the number of 
                            #time steps specified by the variable "steps". 
        
    return data
        

def three_D(L, N, t, steps):
    """
    3-D simulation  
    Assumptions: Particles must move by one cell in each trial;
        they won't stop moving just because they are at the boundary.
        When pairs of particles occupy the same space, they annihilate.
    Boundary conditions for 3-D lattice: 
        If index for any coordinate is 0 or L-1, when the particle would have 
        "crossed over" the boundary, they instead move away from the boundary..
    Inputs: 
        L = length of LxLxL grid
        N = number of initial walkers
        t = time steps for which we will run the simulation
        Steps: number of time steps between which we will make measurements of 
        concentration
    
    Outputs measurements of concentration (of walkers) over time. 
    """
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
        #if there is already a random walker occupying a certain location, 
        #new coordinates for x, y, and z are generated until the walker is given
        #a unique starting position
    
        coord += [[x, y, z]] #new coordinate concatenates into list
    
    coordcopy = np.array(coord)
    #list is converted to numpy array so array operations can be performed
    
    
    for n in range (0, int(t/steps)):
        for j in range (0, steps):
            
            if (coordcopy.size != 0):
                
                for k in range (0, coordcopy[:, 0].size):
                #loops through the coordinates of the random walkers
                #if there are no more walkers left, the loop is not executed
                    
                    chance = (nprnd.random(4))  
                    # generates 4 random floating points numbers from 0 to 1
                    
                    
                    #scenario where particle moves in x, y, and z directions
                    if (chance[3] < (4.0/13)):
                        indices = [0,1,2]
                    #scenario where particle moves in only the x and y directions
                    elif (chance[3] < (6.0/13)): 
                        indices = [0,1]
                    #scenario where particle moves in only the x and z directions   
                    elif (chance[3] < (8.0/13)):
                        indices = [0,2]
                    #scenario where particle moves in only the y and z directions
                    elif (chance[3] < (10.0/13)):
                        indices = [1,2]
                    #scenario where particle moves in only the x direction
                    elif (chance[3] < (11.0/13)):
                        indices = [0]
                    #scenario where particle moves in only the y direction
                    elif (chance[3] < (12.0/13)):
                        indices = [1]
                    #scenario where particle moves in only the z direction
                    else:
                        indices = [2]
                    
                    #implements boundary conditions
                    #loops through the values in the list "indices," so that for 
                    #the coordinate(s) along which the particle can move (determined
                    #by a random number), each coordinate(s) can increment or decrement 
                    #depending on separate random number or its location if it is bordering the boundary
                    #Each coordinate is incremented/decremented independently from the other.
                    #This method allows for scalability when changing the dimensions
                    for m in indices: 
                        if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                            coordcopy[:, m][k] -= 1
                                
                        else: 
                            coordcopy[:, m][k] += 1
            
            #implements the annihilation by updating coordcopy
            #removes pairs of coordinates that are duplicates - ie. pairs of random walkers in the same space 
            coordcopy, count = remove_dup(coordcopy)
         
        data += [count/(L**3)] #records the concentration of walkers after the number of 
                            #time steps specified by the variable "steps". 
        
    return data


def remove_dup(array):
    """
    Function to remove pairs of elements that are duplicates in a numpy array. 
    Implements the annihilation of random walkers. 
    Returns a new numpy array with those elements removed and the number of 
    remaining particles.
    If there are no duplicate elements in the input array, returns the input array.
    If there are no more elements in the array, returns empty array.
    """
    if (array.size > 0): 
        coparray, dupe = np.unique(array, axis = 0, return_counts = True)
        #"dupe" - an array containing the number of occurences for each element
        #in the input array -indicates where the duplicates exist in the input array
        #takes into account the entire set of coordinates for a random walker,
        #not just individual numbers
        #"coparray" returns all unique elements from the input array
        dupe2 = []
            
        #finds indices of elements that had duplicate(s) in the input array 
        #and occurred in an even number 
        for i in range(0, dupe.size):
            if dupe[i] > 1 and ((dupe[i] % 2) == 0):
                dupe2 += [i]
        
        #removes elements in "coparray" that had one or more duplicates in the input array
        #and had an even number of occurences, since two particles are required
        #for annihilation to occur.      
        coparray = np.delete(coparray, dupe2, 0)
        if (len(coparray.shape) == 1): #if the input array is 1D
            size = coparray.size
        else:    #if the input array is more than 1 dimensions
            size = coparray[:,0].size
           
        
    else: #If the input array is empty, indicating that there are no more random walkers.
        size = 0 
        a = []
        coparray = np.array(a)
        
    return coparray, size

#Individual simulations used for further verification.
#These are able to show all the random walkers annihilated, 
    #although there will be a difference in the actual time required between 
    #individual simulations.
#In the Variable explorer, each increment of decrease is 2/(number of cells) 
    #ie. 2/L for 1D, 2/(L**2) for 2D, and 2/(L**3) for 3D
oneDindiv = one_D (40, 10, 5000, 1)
twoDindiv = two_D (10, 30, 300, 1) 
threeDindiv = three_D (4, 30, 100, 1)
    


######## Part 2. Running Multiple Simulations
    
def multisim(func, L, N, t, steps, sims):
    """
    Performs multiple, statistically indipendent simulations
    Inputs: 
        func is the function that runs an individual simulation
        L is the length of one side of the lattice 
        N is the number of initial walkers
        t is the number of time steps for which each simulation will run
        steps: number of time steps between which we will make measurements of 
        concentration
        sims: number of separate simulations
        
    Outputs all of the data from individual simulations in one array 
    and an array containing average values across simulations
    """
    data0 = func(L,N,t,steps)
    for j in range (0, sims):
        data = func(L, N, t,steps)
        data0 = np.vstack([data0, data]) 
    #each row of data0 now contains data from a separate individual simulation
    
    avg = np.mean(data0, axis = 0)
    #array containing averages of values at each time step
    
    return data0, avg

#returns averages across 10 statistically separate simulations for each dimension
oneDrun = multisim(one_D,729, 300, 500, 1, 10)[1] 
twoDrun = multisim(two_D, 27, 300, 180, 1, 10)[1] 
threeDrun = multisim(three_D, 9, 300, 80, 1, 10)[1] 

#plotting averages of [A] across multiple simulations - 1D 
t1 = (np.arange(0,len(oneDrun)))  
plt.plot(t1, oneDrun)
plt.title("1D Lattice: [A] vs Time") 
plt.xlabel("time")
plt.ylabel("[A]")
plt.figure()

#plotting averages of [A] across multiple simulations - 2D 
t2 = (np.arange(0,len(twoDrun)))
plt.plot(t2, twoDrun)
plt.title("2D Lattice: [A] vs Time")
plt.xlabel("time")
plt.ylabel("[A]")
plt.figure()

#plotting averages of [A] across multiple simulations - 3D 
t3 = (np.arange(0,len(threeDrun)))
plt.plot(t3, threeDrun)
plt.title("3D Lattice: [A] vs Time")
plt.xlabel("time")
plt.ylabel("[A]")
plt.figure()




######## Part 3. Plotting with Error Bars
def stdev (indiv, avgs, sims):
    """
    Inputs: 
        indiv - an array where each row contains data for a separate simulation
        avgs - a 1D array containing average values across multiple simulations
        sims - number of statistically independent simulations
        
    Returns a list with the standard deviations of data from 
    separate simulations for each time step.
    """
    sigmas = []
    
    for m in range (0, avgs.size): 
        sumsi = 0
        
        for n in range (0, sims):
            sumsi += ((indiv[n][m]-avgs[m])**2) 
            #sum of the squares of the differences between values from individual simulations and the average
            
        sdev = np.sqrt(sumsi/(sims - 1))
        #sdev is the standard deviation 
        
        sigmas += [sdev]
        
    return sigmas


data1, oneDavg = multisim(three_D, 3, 23, 30, 1, 20)


sigma1 = stdev(data1, oneDavg, 5) #gets standard deviations
t3err = np.arange(0, len(oneDavg)) #creates time array of the proper length

#plots averages of [A] accross multiple simulations, along with error bars
plt.errorbar(t3err, oneDavg, sigma1) 
plt.xlabel("time")
plt.ylabel("[A]")
plt.title("[A] over Time - Error Bars")
plt.figure()




######## Part 4. Graphical Representation of Random Walkers

#The goal in this part is to be able to visualize up to two random walkers. 
#Since we are using a very small number of random walkers, they do not always
#annihilate within the given number of time steps. Sometimes, the code has to be run again 
#or the number of time steps has to be increased, to see the annihilation occur 
#within the given number of time steps.
def two_Dgraphic(L, N, t, steps):
    """
    Almost identical to the function "two_D", but instead of simply deleting 
    elements for which there are duplicates, new coordinates for the random walkers
    are appended to an array. 
    
    Returns all the coordinates that the random walkers would have "traversed."
    """
    coord = []
    
    for i in range (0, N):
        x = rd.randrange(0, L)  #generate random integers - starting positions
        y = rd.randrange(0, L)  #x and y coordinates
        while [x,y] in coord:
            x = rd.randrange(0, L)  
            y = rd.randrange(0, L)
            
        coord += [[x, y]]
    
    coordcopy = np.array(coord)
    coordcopy0 = np.copy(coordcopy)
    
    for n in range (0, int(t/steps)):
        for j in range (0, steps):
            
            if (coordcopy.size != 0):
                
                for k in range (0, coordcopy[:, 0].size):
                    
                    chance = (nprnd.random(3))  # 3 random floats from 0 to 1
                    #print(chance)
                    #scenario where particle moves in both x and y directions
                    if (chance[2] < 0.5):
                        indices = [0, 1]
                    #scenario where particle moves in only the x direction
                    elif (chance[2] < 0.75): 
                        indices = [0]
                    #scenario where particle moves in only the y direction   
                    else:
                        indices = [1]
                        
                    for m in indices: 
                        if ((chance[m] < 0.5) and (coordcopy[:, m][k] > 0)) or (coordcopy[:, m][k] == L - 1): 
                            coordcopy[:, m][k] -= 1
                                
                        else: 
                            coordcopy[:, m][k] += 1
            
            #appends a copy of the new coordinates onto array with old coordinates
            coordcopy0 = np.hstack([coordcopy0, coordcopy])
            
            #implements the annihilation by updating coordcopy
            #removes pairs of coordinates that are duplicates 
            coordcopy, count = remove_dup(coordcopy)
             
            if count == 0:  #If no walkers remain, exits both for loops
                break
        if count == 0:
                break    
        
    return coordcopy0


twoRun = two_Dgraphic (4, 2, 10, 1)


#The following organizes the data collected from the random walks -
#designates a sequence of coordinates to specific corresponding random walkers
index = np.arange(1, int(twoRun[0].size), 2)

coord1 = np.ones([int(twoRun[0].size/2),2]) 
coord2 = np.ones([int(twoRun[0].size/2),2])

row = 0

for i in index:
    coord1[row][1] = twoRun[0][i] #coordinates for one random walker
    coord1[row][0] = twoRun[0][i-1]
    
    coord2[row][1] = twoRun[1][i] #coordinates for the second random walker
    coord2[row][0] = twoRun[1][i-1]
    
    row += 1

#Plots a visual representation of the path of the random walkers
#coordinates of the random walkers are printed in chronological order, so that it 
#is possible to follow the paths of the random walkers. 
print ("Coordinates of Random Walk 1: Blue (in the graphic representation)") 
print(coord1)  
print("Coordinates of Random Walk 2: Red (in the graphic representation)")
print (coord2)
plt.plot(coord1[:,0], coord1[:,1], 'bo-')
plt.plot(coord2[:,0], coord2[:,1], 'ro-')   
plt.figure()



######## Part 5. Performing Linear Fits to Data

def datainv(array):
    """
    For all elements in the input array that are greater than zero,
    returns a list with 1 over those elements.
    """
    data = []
    for i in array:
        if i > 0:
            data += [1/i]
    return data

oneDrun2 = multisim(one_D,729, 650, 500, 1, 10)[1]
twoDrun2 = multisim(two_D, 27, 650, 180, 1, 10)[1]
threeDrun2 = multisim(three_D, 9, 650, 80, 1, 10)[1]

#returns time and 1/concentration values spliced to exclude small concentrations 
x1, y1 = t1[0:700], datainv(oneDrun2)[0:700]
x2, y2 = t2[0:150], datainv(twoDrun2)[0:150]
x3, y3 = t3[0:70], datainv(threeDrun2)[0:70]

#returns slope, y-intercept, and r of linear fit
m1, b1, r1 = sp.linregress(x1, y1)[0:3]
m2, b2, r2 = sp.linregress(x2, y2)[0:3]
m3, b3, r3 = sp.linregress(x3, y3) [0:3]
  

#plotting of 1/[A] and line of best fit - 1D
print("Linear Fit 1: r^2")
print (r1**2) #prints r^2 value
plt.title("1D Lattice: 1/[A] vs Time")
oneDinv = datainv(oneDrun2)
t1a = (np.arange(0,len(oneDinv)))  
plt.plot(t1a[0:700], oneDinv[0:700], 'bo', markersize = 2) #plots data
plt.plot(t1a[0:700], m1*t1a[0:700] + b1, 'r') #plots linear fit
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()

#plotting of 1/[A] and line of best fit - 2D
print("Linear Fit 2: r^2")
print (r2**2) #prints r^2 value
plt.title("2D Lattice: 1/[A] vs Time")
twoDinv = datainv(twoDrun2)
t2a = (np.arange(0,len(twoDinv)))
plt.plot(t2a[0:150], twoDinv[0:150], 'bo', markersize = 2) #plots data
plt.plot(t2a[0:150], m2*t2a[0:150] + b2, 'r') #plots linear fit
plt.xlabel("time")
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()

#plotting of 1/[A] and line of best fit - 3D
print("Linear Fit 3: r^2")
print (r3**2) #prints r^2 value
plt.title("3D Lattice: 1/[A] vs Time")
threeDinv = datainv(threeDrun2)
t3a = (np.arange(0,len(threeDinv)))  
plt.plot(t3a[0:70], threeDinv[0:70], 'bo', markersize = 3) #plots data
plt.plot(t3a[0:70], m3*t3a[0:70] + b3, 'r') #plots linear fit
plt.xlabel("time")
plt.ylabel("1/[A]")
plt.figure()