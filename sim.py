from polymer import Polymer
from particles import Particles
from util import periodic_dist, find_bond_angle, get_positions
from scipy.spatial.distance import pdist, cdist
import matplotlib.pyplot as plt
import numpy as np
from math import *

class Sim:

    def __init__(self, boxDimensions, numMonomers, monomerDiam, particleDiam, particleSpacing,latticeType):
        self.x = boxDimensions[0]
        self.y = boxDimensions[1]
        self.z = boxDimensions[2]

        self.polymer = Polymer(numMonomers, monomerDiam)
        self.particles = Particles(particleDiam, particleSpacing, latticeType, boxDimensions)

        # CONSTANTS FROM PAPER
        self.epsilon = 1
        self.delta = 0.2*monomerDiam
        self.r_eq = 0.8*monomerDiam
        self.r_max = 1.3*monomerDiam
        self.k_f = 100*self.epsilon/(monomerDiam**2)
        self.k_bT = 1

        self.k_theta = 1
        self.epsilon_PN = 1

        if latticeType == 'repulsive':
            self.r_c = 2.5*monomerDiam
        elif latticeType == 'attractive':
            self.r_c = monomerDiam
        else:
            print("Invalid lattice type provided! It must be either attractive or repulsive.")
            exit(-1)


    def initializePolymer(self):
        monomers = []

        boundaries = np.array([self.x, self.y, self.z])

        # We'll start by initializing the first point
        passesChecks = False
        while not passesChecks:
            # Choose a random point in the box
            currentPoint = np.array([np.random.random()*self.x, np.random.random()*self.y, np.random.random()*self.z])

            # Compute the distance between p and each point in q
            distToCurrent = currentPoint - self.particles.lattice

            # Apply the minimum image convention
            distToCurrent = distToCurrent - boundaries * np.round(distToCurrent / boundaries)

            # Need to verify that we're not intersecting with any of the lattice particles
            # Compute the distance between p and each point in q
            distToCurrent = np.linalg.norm(distToCurrent, axis=1)

            # Compute the sum of the radii
            overlap_dist = self.polymer.monomerDiam/2 + self.particles.diameter/2

            # Check if there is any overlap
            overlap = np.any(distToCurrent < overlap_dist)

            if not overlap:
                passesChecks = True

        monomers.append(currentPoint)

        for i in range(self.polymer.numMonomers-1):
            lastPoint = currentPoint
            passesChecks = False

            while not passesChecks:
                # Choose two random angles
                theta = np.random.random()*2*pi
                phi = np.random.random()*2*pi

                # Convert angles to x, y, z, using spherical coordinates
                currentPoint = (lastPoint + np.array([self.polymer.monomerDiam * sin(theta) * cos(phi), 
                                                      self.polymer.monomerDiam * sin(theta) * sin(phi), 
                                                      self.polymer.monomerDiam * cos(theta)])) % boundaries

                # Convert our list of monomers to a numpy array so following operations will work
                monomersSoFar = np.array(monomers)

                distToCurrentForLattice = currentPoint - self.particles.lattice
                distToCurrentForMonomers = currentPoint - monomersSoFar

                # Adjust distances appropriately using periodic boundary conditions
                distToCurrentForLattice = distToCurrentForLattice - boundaries * np.round(distToCurrentForLattice / boundaries)
                distToCurrentForMonomers = distToCurrentForMonomers - boundaries * np.round(distToCurrentForMonomers / boundaries)

                distToCurrentForLattice = np.linalg.norm(distToCurrentForLattice, axis=1)
                distToCurrentForMonomers = np.linalg.norm(distToCurrentForMonomers, axis=1)

                # Calculate distance where there's overlaps
                lattice_overlap_dist = self.polymer.monomerDiam/2 + self.particles.diameter/2
                monomer_overlap_dist = self.polymer.monomerDiam

                # Check if there is any overlap with both the lattice particles and current monomers
                lattice_overlap = np.any(distToCurrentForLattice <= lattice_overlap_dist)
                monomer_overlap = np.any(distToCurrentForMonomers <= monomer_overlap_dist)

                # If it checks out, we can go ahead and add the point
                if not lattice_overlap and not monomer_overlap:
                    passesChecks = True
                
            monomers.append(currentPoint)
        
        self.polymer.monomers = np.array(monomers)
            


    def betweenMonomers(self, monomers):
        E = 0
        boundaries = np.array([self.x, self.y, self.z])

        # Find the pairwise distances between all monomers using pdist w. custom dist function
        pairwiseDists = pdist(monomers, metric=periodic_dist, L=boundaries)

        # Separate bonded vs. non-bonded distances
        bondedInds = get_positions(self.polymer.numMonomers)
        unbondedInds = np.array(list(set(range(self.polymer.numMonomers)) - set(bondedInds)))
        bondedDists = pairwiseDists[bondedInds]
        unbondedDists = pairwiseDists[unbondedInds]

        unbondedOverlaps = unbondedDists[unbondedDists < self.polymer.monomerDiam]
        E += self.epsilon*((self.polymer.monomerDiam/unbondedOverlaps)**12 - (self.polymer.monomerDiam/unbondedOverlaps)**6 + 1).sum()
        
        bondedOverlaps = bondedDists[(bondedDists <= self.r_max) & (bondedDists >= (2*self.r_eq - self.r_max))]
        E += (-self.k_f/2)*(self.r_max - self.r_eq)**2*np.log(1 - ((bondedOverlaps-self.r_eq)/(self.r_max-self.r_eq))**2).sum()
        if len(bondedOverlaps) != len(bondedDists):
            E += inf

        # Generate mx3 array of points that go into all the bond angles
        # Get bond angles, then apply function to get energy there
        groups = np.column_stack((monomers[:-2], monomers[1:-1], monomers[2:]))
        E += ((1/2)*self.k_theta*(pi - np.apply_along_axis(find_bond_angle, 1, groups))**2).sum()


    
        return E
    


    def betweenMonomersandLattice(self, monomers):
        
        boundaries = np.array([self.x, self.y, self.z])
        U_c = -1*self.epsilon_PN*((self.polymer.monomerDiam/self.r_c)**12 - 2*(self.polymer.monomerDiam/self.r_c)**6)

        # Use cdist to get all pairwise distances between lattice and polymer
        dist = cdist(monomers, self.particles.lattice, metric=periodic_dist, L=boundaries)

        overlap_dist = self.polymer.monomerDiam/2 + self.particles.diameter/2

        overlaps = dist[dist < overlap_dist]

        E = (self.epsilon_PN*((self.polymer.monomerDiam/overlaps)**12 
                                - 2*(self.polymer.monomerDiam/overlaps)**6) + U_c).sum()
        
        
        return E


    def run(self, numSteps):
        # Calculate starting energy of the system
        prevE = self.betweenMonomers(self.polymer.monomers) + self.betweenMonomersandLattice(self.polymer.monomers)
        
        # Find the starting center of mass
        startingCenterOfMass = np.mean(self.polymer.monomers, axis=0)

        # Compute the starting radius of gyration
        radiusOfGyration =  np.sqrt(np.sum(np.sum((self.polymer.monomers-startingCenterOfMass)**2, axis=1)) / self.polymer.numMonomers)
    
        startingMonomers = self.polymer.monomers.copy()
        step = 0
        # While we haven't exceeded the center of mass thing, do MC
        while step < numSteps:
            E = 0
            index = np.random.randint(self.polymer.numMonomers)
            currentMonomers = self.polymer.monomers.copy()

            # Move our chosen monomer by dx, dy, dz (which are all somewhere between -delta and +delta)
            currentMonomers[index] = (currentMonomers[index] + (np.array([np.random.random()*2*self.delta, 
                                                                          np.random.random()*2*self.delta, 
                                                                          np.random.random()*2*self.delta]) 
                                                                          -self.delta)) % np.array([self.x, self.y, self.z])
            
            # Get energy for change
            E = self.betweenMonomers(currentMonomers) + self.betweenMonomersandLattice(currentMonomers)

            # print(prevE)
            # print(E)
            # print("-----")

            p = min([1, exp(-1*(E-prevE)/self.k_bT)])

            # print(p)

            accept = np.random.random() <= p
            if accept:
                # print("accept")
                # print("-----")
                self.polymer.monomers = currentMonomers
                prevE = E # Update our previous energy

                # Calculate c of m and r of g again
                centerOfMass = np.mean(self.polymer.monomers, axis=0)
                print(np.linalg.norm(abs(startingCenterOfMass-centerOfMass)))
                radiusOfGyration =  np.sqrt(np.sum(np.sum((self.polymer.monomers-centerOfMass)**2, axis=1)) / self.polymer.numMonomers)
                print(radiusOfGyration)
            # else:
            #     print("reject")
            #     print("-----")
            step += 1

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        polyX = [row[0] for row in self.polymer.monomers]
        polyY = [row[1] for row in self.polymer.monomers]
        polyZ = [row[2] for row in self.polymer.monomers]

        startingX = [row[0] for row in startingMonomers]
        startingY = [row[1] for row in startingMonomers]
        startingZ = [row[2] for row in startingMonomers]


        ax.plot3D(polyX,polyY,polyZ)
        ax.plot3D(startingX, startingY, startingZ)
        plt.show()
                


    def visualize(self, showLattice):
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        polyX = [row[0] for row in self.polymer.monomers]
        polyY = [row[1] for row in self.polymer.monomers]
        polyZ = [row[2] for row in self.polymer.monomers]

        latticeX = [row[0] for row in self.particles.lattice]
        latticeY = [row[1] for row in self.particles.lattice]
        latticeZ = [row[2] for row in self.particles.lattice]

        ax.plot3D(polyX,polyY,polyZ)
        ax.scatter3D(polyX,polyY,polyZ,s=10)
        if showLattice:
            ax.scatter3D(latticeX, latticeY, latticeZ,s=5)
        plt.show()
