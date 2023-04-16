from polymer import Polymer
from particles import Particles
import matplotlib.pyplot as plt
import numpy as np
from math import *

class Sim:

    def __init__(self, boxDimensions, numMonomers, monomerRadius, particleRadius, particleSpacing,latticeType):
        self.x = boxDimensions[0]
        self.y = boxDimensions[1]
        self.z = boxDimensions[2]

        self.polymer = Polymer(numMonomers, monomerRadius)
        self.particles = Particles(particleRadius, particleSpacing, latticeType, boxDimensions)

        # CONSTANTS FROM PAPER
        self.epsilon = 1
        self.delta = 0.2*monomerRadius
        self.r_eq = 0.8*monomerRadius
        self.r_max = 1.3*monomerRadius
        self.k_f = 100*self.epsilon/(monomerRadius**2)
        self.k_bT = 1

        self.k_theta = 1
        self.epsilon_PN = 1

        if latticeType == 'repulsive':
            self.r_c = 2.5*monomerRadius
        elif latticeType == 'attractive':
            self.r_c = monomerRadius
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
            overlap_dist = self.polymer.monomerRadius + self.particles.radius

            # Check if there is any overlap
            overlap = np.any(distToCurrent < overlap_dist)

            if not overlap:
                passesChecks = True

        monomers.append(currentPoint)

        for i in range(self.polymer.numMonomers-1):
            lastPoint = currentPoint
            passesChecks = False
            while not passesChecks:

                theta = np.random.random()*2*pi
                phi = np.random.random()*2*pi

                currentPoint = (lastPoint + np.array([self.polymer.monomerRadius* 2 * sin(theta) * cos(phi), 
                                                      self.polymer.monomerRadius * 2 * sin(theta) * sin(phi), 
                                                      self.polymer.monomerRadius * 2 * cos(theta)])) % boundaries

                monomersSoFar = np.array(monomers)
                

                distToCurrentForLattice = currentPoint - self.particles.lattice
                distToCurrentForMonomers = currentPoint - monomersSoFar

                distToCurrentForLattice = distToCurrentForLattice - boundaries * np.round(distToCurrentForLattice / boundaries)
                distToCurrentForMonomers = distToCurrentForMonomers - boundaries * np.round(distToCurrentForMonomers / boundaries)

                distToCurrentForLattice = np.linalg.norm(distToCurrentForLattice, axis=1)
                distToCurrentForMonomers = np.linalg.norm(distToCurrentForMonomers, axis=1)

                lattice_overlap_dist = self.polymer.monomerRadius + self.particles.radius
                monomer_overlap_dist = self.polymer.monomerRadius*2

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

        for i in range(len(monomers)):
            # Find the energy associated with non-bonded monomer interactions
            bondedMonomers = []
            nonBondedMonomers = []

            # Find the energy associated with bonded monomer interactions

            r = 0

            if r <= self.r_max and r >= (2*self.r_eq - self.r_max):
                E += 1
            else:
                E += inf

        

        return E
    
    def betweenMonomersandLattice(self, monomers):
        # Find the energy associated with monomers close to lattice points
        for monomer in monomers:
            print("yes")
        return 0
    
    def run(self):
        # Calculate starting energy
        E = self.betweenMonomers(self.polymer.monomers) + self.betweenMonomersandLattice(self.polymer.monomers)
        
        centerOfMass = np.mean(self.polymer.monomers, axis=0)

        # Compute the radius of gyration
        radiusOfGyration =  np.sqrt(np.sum(np.sum((self.polymer.monomers-centerOfMass)**2, axis=1)) / self.polymer.numMonomers)
        

        # Choose a random monomer to displace
        prevE = E
        # While we haven't exceeded the center of mass thing, do MC
        while True:
            index = np.random.randint(self.polymer.numMonomers)
            currentMonomers = self.polymer.monomers

            # Move our chosen monomer by dx, dy, dz (which are all somewhere between -delta and +delta)
            currentMonomers[index] = (currentMonomers[index] + (np.array([np.random.random()*2*self.delta, 
                                                                          np.random.random()*2*self.delta, 
                                                                          np.random.random()*2*self.delta]) 
                                                                          -self.delta)) % np.array([self.x, self.y, self.z])
            
            # Get energy for change
            E = self.betweenMonomers(currentMonomers) + self.betweenMonomersandLattice(currentMonomers)

            p = min(1, exp(-1*(E-prevE)/self.k_bT))

            accept = np.random.random() <= p
            if accept:
                self.polymer.monomers = currentMonomers
                preV = E

            # Calculate c of m and r of g again

    

    
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