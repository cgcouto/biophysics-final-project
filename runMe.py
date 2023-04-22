from sim import Sim

boxDimensions = [100,100,100]
numMonomers = 20
monomerDiam = 2
particleDiam = 2
particleSpacing = 10
latticeType = 'repulsive'

sim = Sim(boxDimensions, numMonomers, monomerDiam, particleDiam, particleSpacing, latticeType)

sim.initializePolymer()

# sim.visualize(False)

sim.run(100)

# sim.visualize(False)


