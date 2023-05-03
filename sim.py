from scipy.spatial.distance import pdist
from util import find_bond_angle, get_positions
import matplotlib.pyplot as plt
import numpy as np
from math import *
import statistics
import scipy.stats as stats

# Find the starting polymer chain within a given vesicle
# numMonomers (int) : an integer number of monomers to place
# monomerDiam (float) : the diameter of each monomer
# vesicleDiam (float) : the diameter of the vesicle
def initializePolymer(numMonomers, monomerDiam, vesicleDiam):
    if vesicleDiam < monomerDiam:
        print("ERROR: The diameter of the vesicle cannot be smaller than the diameter of a monomer.")
        exit(-1)

    # We'll start by initializing the first point
    passesChecks = False
    while not passesChecks:

        # Choose a random point in a box with side lengths equal to vesicleDiam
        currentPoint = np.random.rand(1,3)*vesicleDiam - vesicleDiam/2
        distFromCenter = np.linalg.norm(currentPoint)

        # Accept if it's within the vesicle
        if distFromCenter <= vesicleDiam/2:
            passesChecks = True

    monomers = currentPoint
    print("monomer 1 placed")

    trueStep = 1
    while len(monomers) != numMonomers:
        lastPoint = currentPoint.copy()
        passesChecks = False

        step = 0

        while not passesChecks:
            # Choose two random angles
            theta = np.random.random()*2*pi
            phi = np.random.random()*2*pi

            # Convert angles to x, y, z, using spherical coordinates
            currentPoint = (lastPoint + np.array([monomerDiam * sin(theta) * cos(phi), 
                                                    monomerDiam * sin(theta) * sin(phi), 
                                                    monomerDiam * cos(theta)]))

            distFromCenter = np.linalg.norm(currentPoint)
            distToOtherMonomers = np.linalg.norm(monomers - currentPoint, axis=1)

            # Check if there is any overlap with current monomers and that we're within the vesicle
            monomerOverlap = np.sum(distToOtherMonomers <= monomerDiam) > 0
            insideVesicle = distFromCenter <= vesicleDiam/2

            # If it checks out, we can go ahead and add the point
            if not monomerOverlap and insideVesicle:
                passesChecks = True
                trueStep += 1

            step += 1

            # Need this because the monomer placing process can get stuck...
            # If it takes too many steps we restart from the first placed monomer
            if step > 250000:
                print("restart")
                monomers = np.delete(monomers, np.s_[1:], axis=0) # Restart to just the first point again
                lastPoint = monomers.copy()
                step = 0
                trueStep = 1
            
        monomers = np.append(monomers, currentPoint.reshape(1, 3), axis=0)
        print("monomer " + str(trueStep) + " placed")
        
    return monomers
    
        
# Calculate the potential interactions that occur between monomers
# (this included bonded dists, bonded angles, and unbonded dists)
# monomers (nx3 numpy array) : monomer positions
# numMonomers (int) : number of monomers in the chain
# monomerDiam (float) : the diameter of the monomers
# k_theta (int) : chain stiffness coefficient
def betweenMonomers(monomers, numMonomers, monomerDiam, k_theta):
    # Starting energy should be zero!
    E = 0

    # Constants from Peng et al. that I didn't change at all
    epsilon = 1
    r_eq = 0.8*monomerDiam
    r_max = 1.3*monomerDiam
    k_f = 100*epsilon/(monomerDiam**2)

    # Find the pairwise distances between all monomers using pdist
    pairwiseDists = pdist(monomers)

    # Separate bonded vs. non-bonded distances
    bondedInds = get_positions(numMonomers)
    unbondedInds = np.array(list(set(range(len(pairwiseDists))) - set(bondedInds)))
    bondedDists = pairwiseDists[bondedInds]
    unbondedDists = pairwiseDists[unbondedInds]

    # Find which of the unbonded distances are too small, add to energy accordingly
    unbondedOverlaps = unbondedDists[unbondedDists < monomerDiam]
    E += epsilon*((monomerDiam/unbondedOverlaps)**12 - (monomerDiam/unbondedOverlaps)**6 + 1).sum()
    
    # Find which of the bonded distances leave the sweet spot, add to energy accordingly
    bondedOverlaps = bondedDists[(bondedDists <= r_max) & (bondedDists >= (2*r_eq - r_max))]
    E += (-1*k_f/2)*(r_max - r_eq)**2*np.log(1 - ((bondedOverlaps-r_eq)/(r_max-r_eq))**2).sum()
    
    if len(bondedOverlaps) != len(bondedDists):
        E += inf

    # Generate mx3 array of points that go into all the bond angles
    # Get bond angles, then apply function to get energy there
    groups = np.column_stack((monomers[:-2], monomers[1:-1], monomers[2:]))
    E += ((1/2)*k_theta*(pi - np.apply_along_axis(find_bond_angle, 1, groups))**2).sum()

    return E

# We only check the moved monomer each time step as that's the only thing that COULD be outside!
# changedMonomer (1x3 numpy array) : the x-y-z position of the monomer that's been moved
# vesicleDiam (float) : the diameter of the simulated vesicle
def constrainedInVesicle(changedMonomer, vesicleDiam):

    # Find the distance from the monomer to the center of the vesicle
    distFromCenter = np.linalg.norm(changedMonomer)

    # If we're within we're good, otherwise reject state with inf energy
    if distFromCenter <= vesicleDiam/2:
        E = 0
    else:
        E = inf
    return E


# Runs the Monte Carlo simulation!
# numSteps (int) : the number of time steps to run for
# monomers (nx3 numpy array) : the starting polymer positions
# numMonomers (int) : number of monomers in the chain
# monomerDiam (float) : the diameter of the monomers
# vesicleDiam (float) : the diameter of the vesicle
# k_theta (int) : chain stiffness coefficient
def run(numSteps, monomers, numMonomers, monomerDiam, vesicleDiam, k_theta):
    
    # Additional constants from Peng et al. that I didn't adjust
    delta = 0.2*monomerDiam
    k_bT = 1

    # Start the simulation data with the initial configuration
    frameData = []
    frameData.append(monomers)

    # Initialized boolean timestep array
    movedMonomers = np.array([False]*numMonomers)

    # Calculate starting energy of the system
    # no constrainedInVesicle as initializePolymer made sure everything's in there
    prevE = betweenMonomers(monomers, numMonomers, monomerDiam, k_theta)

    step = 0

    while step < numSteps:
 
        index = np.random.randint(numMonomers) # Choose random monomer to displace
        currentMonomers = monomers.copy()

        # Move our chosen monomer by dx, dy, dz (which are all somewhere between -delta and +delta)
        currentMonomers[index] = currentMonomers[index] + (np.array([np.random.random()*2*delta, 
                                                                        np.random.random()*2*delta, 
                                                                        np.random.random()*2*delta]) 
                                                                        -delta)
        
        # Get energy for change
        E = betweenMonomers(currentMonomers, numMonomers, monomerDiam, k_theta) + constrainedInVesicle(currentMonomers[index], vesicleDiam)

        # Find probability for accepting new configuration
        p = min([1, exp(-1*(E-prevE)/k_bT)])

        # The moment of truth...
        accept = np.random.random() <= p
        if accept:
            # Update our energy, list of monomers, bool array for timesteps
            monomers = currentMonomers
            prevE = E 
            movedMonomers[index] = True

            # If all the monomers have been moved at least once, increment step count
            # and reset boolean array
            if np.all(movedMonomers):
                step += 1
                print("completed step " + str(step))
                frameData.append(monomers)
                movedMonomers = np.array([False]*numMonomers)

    return frameData
            
# Plot a list of protein configurations in 3D space
# monomersList (list of nx3 numpy arrays) : The protein configurations you wish to plot
def visualize(monomersList):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # Go through the list and plot everything
    for points in monomersList:

        polyX = [row[0] for row in points]
        polyY = [row[1] for row in points]
        polyZ = [row[2] for row in points]

        ax.plot3D(polyX,polyY,polyZ)
        ax.scatter3D(polyX,polyY,polyZ,s=10)
    plt.show()

# Find all three size metrics on a given set of monomers
# monomers (nx3 numpy array) : monomer x-y-z positions
def getSize(monomers):
    endToEndDist = np.linalg.norm(monomers[0]-monomers[-1])
    maxMonomerDist = max(pdist(monomers))
    radiusOfGyration = np.sqrt(np.mean(np.linalg.norm(monomers - np.mean(monomers), axis=1)**2))

    return [endToEndDist, maxMonomerDist, radiusOfGyration*2]

# Computes and plots the protein size as the number of monomers is varied
# numMonomersList (list of ints) : list of varying numbers of monomers in the chain
# monomerDiam (float) : the diameter of the monomers
# vesicleDiam (float) : the diameter of the vesicle
# numTrials (int) : number of trials to run for each numMonomers combo
# k_theta (int) : chain stiffness coefficient
def size_per_num_monomers(numMonomersList, monomerDiam, vesicleDiam, numTrials, k_theta):
    endToEnd = []
    largestDist = []
    diamOfGyration = []
    for count in numMonomersList:
        results = []
        for trial in range(numTrials):
            startingMonomers = []
            startingMonomers = initializePolymer(count, monomerDiam, vesicleDiam)

            endingMonomers = run(100, startingMonomers, count, monomerDiam, vesicleDiam, k_theta)

            results.append(getSize(endingMonomers))

        endToEnd.append([statistics.mean([sizes[0] for sizes in results]), stats.sem([sizes[0] for sizes in results])])
        largestDist.append([statistics.mean([sizes[1] for sizes in results]), stats.sem([sizes[1] for sizes in results])])
        diamOfGyration.append([statistics.mean([sizes[2] for sizes in results]), stats.sem([sizes[2] for sizes in results])])

    print(numMonomersList)
    print(endToEnd)
    print(largestDist)
    print(diamOfGyration)
    plt.errorbar(numMonomersList, [data[0] for data in endToEnd], yerr= [data[1] for data in endToEnd], label="end-to-end distance")
    plt.errorbar(numMonomersList, [data[0] for data in largestDist], yerr= [data[1] for data in largestDist], label="largest distance between two monomers")
    plt.errorbar(numMonomersList, [data[0] for data in diamOfGyration], yerr= [data[1] for data in diamOfGyration], label="diameter of gyration")
    plt.legend(loc ='upper left')
    plt.title("Polymer size vs. number of monomers, vesicle of diameter {}, chain stiffness emphasized".format(vesicleDiam))
    plt.xlabel("Number of monomers")
    plt.ylabel("Polymer size")

    plt.show()

# Computes and plots the persistence length vs. the mean square end-to-end distance
# numMonomers (ints) : number of  monomers in the chain
# monomerDiam (float) : the diameter of the monomers
# vesicleDiam (float) : the diameter of the vesicle
# numTrials (int) : number of trials to run for each numMonomers combo
# k_thetas (list of ints) : list of varying chain stiffness coefficients
def mean_squared_end_to_end(numMonomers, monomerDiam, vesicleDiam, numTrials, k_thetas):
    R2_sems = []
    R2_temp = []
    for k_theta in k_thetas:
        results = []
        for i in range(numTrials):
            print(k_theta)
            startingMonomers = initializePolymer(numMonomers, monomerDiam, vesicleDiam)
            sim = run(100, startingMonomers, numMonomers, monomerDiam, vesicleDiam, k_theta)
            results.append(np.linalg.norm(sim[-1][-1] - sim[-1][0])**2)
        R2_temp.append(np.mean(results))
        R2_sems.append(stats.sem(results))

    L = numMonomers-1

    R2_actual = np.array(R2_temp)
    Lp_actual = np.array(k_thetas)*0.82*monomerDiam

    Lp_fit = np.linspace(1, 100, 100)*0.82*monomerDiam
    R2_fit = 2*(Lp_fit**2)*((L/Lp_fit)-1+np.exp(-L/Lp_fit)) # Eq. 11 in Peng et al.

    print(Lp_actual)
    print(R2_actual)
    print(R2_sems)

    plt.plot(Lp_fit, R2_fit, label="estimated curve")
    plt.errorbar(Lp_actual, R2_actual, yerr= R2_sems, label="simulated results")
    plt.legend(loc ='upper left')

    plt.xlabel("<Lp>")
    plt.ylabel("<R2>")

    plt.title("Persistence length vs. mean squared end-to-end distance, polymer constrained by vesicle")

    plt.show()