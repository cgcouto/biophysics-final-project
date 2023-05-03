from sim import *

# This is the script where I ran different parts of my code base!

# Set important parameters
numMonomers = 20
monomerDiam = 1
vesicleDiam = 5

simulationSteps = 100
k_theta = 0

## Code for running a single simulation
startingMonomers = initializePolymer(numMonomers, monomerDiam, vesicleDiam)
frameData = run(simulationSteps, startingMonomers, numMonomers, monomerDiam, vesicleDiam, k_theta)

# Visualize the initial and final chains if you wish
visualize([frameData[0], frameData[-1]])

## Code for finding protein size versus the number of monomers
# numMonomersRange = [3,5,10,20,30,40]
# numTrials = 10

# size_per_num_monomers(numMonomersRange, monomerDiam, vesicleDiam, numTrials, k_theta)

## Code for finding persistence length vs. mean square end-to-end distance

# k_thetaRange = [3, 5, 7, 10, 20, 40, 60, 80, 100]
# numTrials = 10

# mean_squared_end_to_end(numMonomers, monomerDiam, vesicleDiam, numTrials, k_thetaRange)
