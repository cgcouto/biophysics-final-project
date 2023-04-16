from sim import Sim

boxDimensions = [100,100,100]
numMonomers = 20
monomerRadius = 1
particleRadius = 2
particleSpacing = 10
latticeType = 'repulsive'

sim = Sim(boxDimensions, numMonomers, monomerRadius, particleRadius, particleSpacing, latticeType)

sim.initializePolymer()

sim.visualize(False)
