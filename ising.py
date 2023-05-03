#!/user/bin/python
import numpy as np
from matplotlib.pylab import *
import time

# Define Parameters
N = 50 # lattice length
T = 0.1 # Temperature
J = 2. # Interaction Energy
h = 0.0 # External field, dramatically increases calculation time if not 0
steps = -1 # number of total steps, inf if negative
update_ix = 1000 # iterations before updating figure

# Initialize spin configuration
L = N**2 # total number of sites on lattice
B = 1./T # inverse temperature

# initialize spin configuration as random configuration w/ value -1 or +1
spin_config = np.random.randint(low=0,high=2,size = (N,N))*2 -1

ion() # set interacting on for matplotlib, neccessary for updating figure
fig = figure('Ising Model') # Create figure
img = imshow(spin_config,cmap='Greys_r',interpolation = 'none') # create image plot
ix = 0 # Set index (iterations) to zero
while (ix < steps) or steps < 0:
    # Select spin with selection probability = 1/L
    row = np.random.randint(N)
    column = np.random.randint(N)

    # Select spins for Hamiltonian interaction
    sij = spin_config[row,column]
    s1 = spin_config[(row+1)%N,column]
    s2 = spin_config[(row-1)%N,column]
    s3 = spin_config[row,(column+1)%N]
    s4 = spin_config[row,(column-1)%N]
    
    # Define Hamiltonian of spin configuration
    H1 = -1*J*sij*(s1+s2+s3+s4)
    
    # Define Hamiltonian with spin flipped
    H2 = J*sij*(s1+s2+s3+s4)  

    # Apply external magnetic field if present -> dramatically slows algorithm
    if h != 0:
        H1 += h*spin_config.sum()
        H2 += h*spin_config.sum() - 2*h*sij

    # Determine Acceptance
    if (H2-H1) <= 0: # lower energy, Accept
        spin_config[row,column] = -1*spin_config[row,column]
    else: # if higher energy, decide based on probability
        probability = np.exp(-1.*B*(H2-H1))
        if np.random.random() < probability:
            spin_config[row,column] = -1*spin_config[row,column]

    # update plot
    if (ix % update_ix) == 0:
        print('steps:', ix) # print progress
        img.set_data(spin_config)
        tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
        tick_params(axis='y',which='both',right='off',left='off',labelleft='off')
        fig.canvas.draw() # Update figure
        fig.canvas.flush_events() # Flush events to update figure

    # increment index (iteration)
    ix += 1
