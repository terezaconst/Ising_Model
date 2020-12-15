from IsingModel import Ising
import matplotlib
matplotlib.use("TkAgg")  # need to set the TkAgg backend explicitly otherwise it introduced a low-level error
from matplotlib import pyplot as plt
matplotlib.rcParams.update({'font.size': 14}) # improve visual presentation for some plots
import scipy as sc
from  scipy import ndimage
from pylab import *
from scipy.optimize import curve_fit
import time
start = time.time()

# Illustrate 2d ferromagnet array
def showIllustration(sample, EquilibrationTime):
    sample.takeTimeStep(timesteps=EquilibrationTime)
    plt.imshow(sample.config)
    plt.xlabel('x direction')
    plt.ylabel('y direction')
    plt.show()


# energy calculations
def heatCapacity(array, temps, time):
    # returns a list containing the heat capacity values for each given temperature.
    # empty arrays to be filled with values from each sample using for loop
    sDeviation = np.zeros((len(array), len(temps)))

    for arrayIndex, sample in enumerate(array):
        for tempIndex, temperature in enumerate(temps):
            sample.temperatureValue(temperature)
            energyRecord = np.zeros(time)
            sample.takeTimeStep(timesteps=1500) # 1500 timesteps to let sample reach equilibrium
            for t in range(time): # energy for each timestep
                energyRecord[t] = sample.energyValue()
                sample.takeTimeStep()
            sDeviation[arrayIndex][tempIndex] = np.std(energyRecord) #standard deviation of energy

    heatCapacities = np.average((sDeviation / temps) ** 2, axis=0)
    return heatCapacities


# calculate critical temperature
def calcCriticalTemp(N, noOfMeasurements, noOfSamples, measurementTime, finiteSizeTest = False, repeats=400):

    minTemp, maxTemp = 1, 5
    samples = [Ising(N, maxTemp) for i in range(noOfSamples)] # samples from ferromagnet class
    temps = np.flipud(np.linspace(minTemp, maxTemp, noOfMeasurements)) # range of temperatures 

    # heatCapacity used to find the average heat capacity of a given number of samples of the same size over many temps
    heatCapacities = heatCapacity(samples, temps, measurementTime)
    # measurements convolved with a Gaussian to smooth out
    convolved = sc.ndimage.gaussian_filter(heatCapacities, sigma=3.2)
    criticalTemp = temps[np.argmax(convolved)]     # the heat capacity peaks at the Curie temperature

    if finiteSizeTest: # produce average value of Tc for a given value of N
        Tc = np.zeros(repeats)
        for i in range(0, repeats):
            Tc[i] = criticalTemp
        Tc_ave = np.average(Tc)
        print('Average T_C = ' +str(Tc_ave))
        return Tc_ave


    else:  # plot C against temps
       plt.figure()
       plt.plot(temps, heatCapacities) 
       plt.plot(temps, convolved)
       plt.xlabel('Temperature / J/k')
       plt.ylabel('Heat Capacity / k')

       print(str(criticalTemp))
       return criticalTemp

# magnetisation tests and plots
def magnetisationTest(sample, sometime, N):
    sample.takeTimeStep(timesteps=1500)  # 1500 timesteps to let sample reach equilibrium
    M = np.zeros(sometime)

    for tau in range(sometime): # M for each timestep
        M[tau] = sample.calcMagnetisation() 
        sample.takeTimeStep()

    print("M[tau] array is  " +str(M))
    plt.plot(np.arange(0, len(M), 1), M/N**2)
    plt.title('Magnetisation against Time')
    plt.xlabel('Number of Iterations')
    plt.ylabel('Magnetisation Per Site')
    plt.show


# calculate and plot autocorrelation with exponential fit
def autoCorrelation(sample, longTime):
    # compute empirical autocovariance with lag tau averaged over time longTime

    sample.takeTimeStep(timesteps=1500)  # 1500 timesteps to let sample reach equilibrium
    M = np.zeros(longTime)
    for tau in range(longTime): # M for each timestep
        M[tau] = sample.calcMagnetisation()
        sample.takeTimeStep()

    M_ave = np.average(M)  # time - average
    M = (M - M_ave)

    autocorrelation = np.correlate(M, M, mode='full')
    autocorrelationArray = autocorrelation[int(len(autocorrelation)/2):]
    autocorrelationArrayforFit = np.absolute(autocorrelationArray/autocorrelationArray[0]) # normalise such that max autocorrelation is 1
    x = np.arange(0, len(autocorrelationArrayforFit), 1)

    # apply exponential fit
    def exponential(x, a, b):
        return a * np.exp(-b * x)

    popt, pcov = curve_fit(exponential, x, autocorrelationArrayforFit)  # array, 2d array
    yy = exponential(x, *popt)
    plt.plot(x, autocorrelationArrayforFit)
    plt.plot(x, yy)
    plt.title('Magnetisation Autocorrelation against Time for Temperature = ' + str(T) + ' J/k')
    plt.ylim(0, 1)
    plt.xlabel('Time / Number of Time Steps ')
    plt.ylabel('Magnetisation Autocorrelation')
    plt.show()

    # prints tau_e value 1/b from exponential a * np.exp(-b * x)
    print('tau_e is ' + str(1/popt[1])) # units converted to time steps by taking reciprocal

#finite size testing
def testFiniteSizeScaling(N, Tinf, a, nu):
    # equation as defined in course handout
    return (Tinf + a * np.power(1. / np.array(N), (1 / nu)))

# plot domain size
def testDomain(sample, testTime):
    sample.takeTimeStep(timesteps=1500)
    showIllustration(sample, EquilibrationTime=1500) # to verify correctness of function
    size = np.zeros(testTime)
    for t in range(testTime):
        size[t] = sample.domainSize()
        sample.takeTimeStep()
    domainSize = np.average(size)
    print( str(domainSize))
    return domainSize

# test hysteretic behaviour
def testHysteresis(sample, rangeH, timePeriod, noOfCycles):

   # empty array to be filled with  magnetisations of each cycle
   results = np.zeros((noOfCycles, len(rangeH)))

   for i in range(noOfCycles):
       for j, H in enumerate(rangeH):
           sample.couplingWithField(H)
           sample.takeTimeStep(timesteps=1500) # timesteps to allow system to reach equilibrium
           magnet = 0
           for t in range(timePeriod):
               magnet += sample.calcMagnetisation() # sum over time
               sample.takeTimeStep()
           mag[i, j] = magnet / timePeriod # magnetisation

   magnetisations = np.average(mag, 0) 
   errors = np.std(results, 0)
   return magnetisations, errors



if __name__ == '__main__':

# plot autocorrelation against time
  longTime = 1000
  N = 30
  temp = [1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 ]
  for T in temp:
      magnet = Ising(30, T)  # (N, temp)
      autoCorrelation(magnet, longTime)

          
# tau_e values obtained from running autocorrelation function for each T above
   tau_e = [2.41, 2.65, 20.74, 39.54, 105.18, 69.86, 21.69, 10.75, 2.32]
   temps = [1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 ]
   plt.plot(temps, tau_e)
   plt.xlabel('Temperature / J/k')
   plt.ylabel('Time lag tau_e')
   plt.title('Time Lag tau_e against Temperature')
   plt.show()


# Creating images of the ferromagnet at equilibrium - performed for coupling = 0, 0.75
	sample2 = Ising(80, 1, coupling = 0.75)
	showIllustration(sample2, EquilibrationTime=2000)
	
	sample25 = Ising(80, 2, coupling = 0.75)
	showIllustration(sample25, EquilibrationTime=2000)
	
	sample3 = Ising(80, 2.5, coupling = 0.75)
	showIllustration(sample3, EquilibrationTime=2000)
	
	sample4 = Ising(80, 4, coupling = 0.75)
	showIllustration(sample4, EquilibrationTime=2000)

    
# plots of magnetisation against time to verify correctness of model
	N = 40
	sometime = 10000 
	temperatures = [1, 2.3, 3] #below, near, above Tc
	for T in temperatures:
	    magnetTest = Ising(N, T, coupling = 0) # completed for coupling = 0 and coupling = 0.75
	    magnetisationTest(magnetTest, sometime, N)
	    plt.show()


# plot heat capacity against temperature and prints critical temperature
	N = 40
	calcCriticalTemp(N, 400, 1, 400, finiteSizeTest = False) #(N,number of measurements ,number of samples, averaging time)
	plt.title('Heat Capacity against Temperature for a square lattice of N = ' + str(N))
	plt.xlim(1, 5)
	plt.show()

# plots of Tc against N and finite size scaling
# values of critical temperatures obtained from running calcCriticalTemp for sizes N
    critTemps = [2.398, 2.3734, 2  3634, 2.3433, 2.3310, 2.3132, 2.3190, 2.3092, 2.2996, 2.2976, 2.2967]
    valueN = [8, 10, 11, 12, 14, 16, 18, 20, 22, 28, 32]

    parameters, errors = sc.optimize.curve_fit(testFiniteSizeScaling, valueN, critTemps) #array, 2d array
    plt.errorbar(valueN, critTemps, np.repeat(4/400, len(valueN)))  # the error = temp range / noOfMeasurements
    plt.plot(valueN, testFiniteSizeScaling(valueN, parameters[0], parameters[1], parameters[2]))
    plt.xlabel('Lattice Width/ N')
    plt.ylabel('Critical Temperature / J/k')
    print(parameters) #determine Tc(inf),a and nu
    print(np.diag(errors)) # diagonals of 2d covariance array provide variance
    plt.title('Critical Temperature Against Lattice Width')
    plt.show()

        
# prints average value of T_c to insert in critTemps array above
	N =[8, 10, 11, 12, 14, 16, 18, 20, 22, 28, 32]
    for Nval in N:
        calcCriticalTemp(Nval, 400, 1, 2000, finiteSizeTest=True)

        
# hysteresis testing - plots for 2 temperatures
	 N = 30
	 temps = [0.2, 3] #below Tc, above Tc
	 numberOfCycles = 4
	 time = 1000
	 Hmin, Hmax = -3.5, 3.5
	 field = np.arange(Hmin,  Hmax, 0.1) # H range with 0.1 increments
	 field = np.append(field, np.flipud(field)) # min to max, and max to min so that it cycles round
	 for T in temps:
	    sample = Ising(N, T)
	    magnetisations, errorBars = testHysteresis(sample, field, time, numberOfCycles)
	    plot = plt.figure()
	    sitemagnetisations = magnetisations / N**2
	    plt.errorbar(field, sitemagnetisations)
	    plt.title('Average Magnetisation per Site against Magnetic Field Strength')
	    plt.xlabel('Magnetic Field Strength/ J/'\u03BC')
	    plt.ylabel('Magnetisation Per Site')
	 plt.show()
	
                   
# domain size testing plots
	N = 30
	field = [0, 0.75]
	temperature = 4 # varied between 0.1 - 4.5 to produce domain size vs temperature plot
	for coupling in field:
	    state = Ising(N, temperature, coupling)
	    domainSize = testDomain(state, testTime = 1000) #testTime must be greater than tau_e
	    domainProportion = domainSize / N**2
	    print('Size of domain as a fraction of the total number of spins is = ' +str(domainProportion))


# plot domain size with temperature
	temperatures = [0.1, 1, 1.5,  2, 2.2, 2.3, 2.4, 2.6, 3, 4, 4.5]
	domainSize = [1, 1, 0.9847, 0.9151, 0.6834, 0.3091, 0.3159, 0.2872, 0.2163, 0.1754, 0.1731] #as a fraction of the total ferromagnet area
	domainSizeField = [1, 1, 0.9960, 0.9726, 0.9594, 0.9521, 0.9404, 0.9060, 0.8380, 0.5640, 0.4670]
	plt.plot(temperatures, domainSize)
	plt.plot(temperatures, domainSizeField)
	plt.xlabel('Temperature / J/k')
	plt.ylabel('Domain Size as Proportion of Total Number of Spins')
	plt.title('Domain Size against Temperature')
	plt.show()

# measure CPU time for each process
end = time.time()
print( 'CPU run time = ' +str(end - start) + ' seconds')
