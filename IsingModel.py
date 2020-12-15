
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # need to set the TkAgg backend explicitly otherwise it introduced a low-level error
from matplotlib import pyplot as plt


# Ising Model Simulation
class Ising():

    # __init__ indicates that the directory should be regarded as a package
    def __init__(self, N, temp, coupling = 0): # coupling is zero unless stated otherwise

        # declarations NB "self" to access the attributes and methods of the class
        self.k = 1  # setting k and J to 1 such that the
        self.J = 1  # temperature is in units of J/k
        self.N = N
        self.temp = temp
        self.energy = 0
        self.magnetisation = 0
        self.coupling = coupling


        # forms the random initial random spin N x N configuration
        self.config = 2 * np.random.randint(2, size=(N, N)) - 1  # -1 to correct spins to (-1,1)


    # Calculate energy and magnetisation of starting configuration.
        for i in range(0, self.N):
            for j in range(0, self.N):
                energyCalc = 0
                spin = self.config[i, j]
                # Periodic boundary conditions used such that every spin has 4 neighbours
                neighbours = self.config[(i + 1) % N, j] + self.config[i, (j + 1) % N] + self.config[(i - 1) % N, j] + \
                             self.config[i, (j - 1) % N]
                # spin-spin interactions
                energyCalc -= spin * neighbours * self.J

                energyCalc -= self.coupling * spin # adding extra H term
                self.magnetisation += spin
                self.energy += energyCalc


    def takeTimeStep(self, timesteps=1):
        for t in range(timesteps):
            for i in range(self.N):
                for j in range(self.N):
                    a = np.random.randint(0, self.N) # only used for magnetisation-related tests
                    b = np.random.randint(0, self.N)
                    # declaration of variables
                    J = self.J
                    N = self.N
                    spin = self.config[i, j]
                    # Periodic boundary conditions used such that every spin has 4 neighbours
                    neighbours = self.config[(i + 1) % N, j] + self.config[i, (j + 1) % N] + self.config[
                        (i - 1) % N, j] + self.config[i, (j - 1) % N] # a/b for magnetisation, i/j for energy
                    energyFlip = J * spin * neighbours
                    energyFlip += self.coupling * spin # adding extra H term

                    if energyFlip < 0:
                        self.energy += 2 * energyFlip  # The change in energy is -2E
                        spin *= -1  # flip spin
                        self.magnetisation += 2 * spin  # The change in magnetisation is +2spin
                        self.config[i,j] = spin
                    elif np.exp(- 2 * energyFlip / (self.k * self.temp)) > np.random.rand():
                        self.energy += 2 * energyFlip
                        spin *= -1 # flip spin
                        self.magnetisation += 2 * spin
                        self.config[i, j] = spin

	# Calculate domain size 
    def domainSize(self):
        self.domain = 0
        N = self.N
        for i in range(self.N):
            for j in range(self.N):
                neighb = self.config[(i + 1) % N, j] + self.config[i, (j + 1) % N] + self.config[(i - 1) % N, j] + \
                         self.config[i, (j - 1) % N]
                if neighb > 0 and self.config[i, j] == self.config[0, 0]:
                    self.domain += 1
        return self.domain




    # Have individual callable functions for each
    def calcMagnetisation(self):
        return self.magnetisation

    def couplingWithField(self, newCoupling):
        self.coupling = newCoupling

    def energyValue(self):
        return self.energy

    def temperatureValue(self, newTemp):
        self.temp = newTemp

