import numpy as np

class Particles:

    def __init__(self, radius, spacing, type, boxDimensions):
        self.radius = radius
        self.spacing = spacing
        self.type = type

        if radius <= spacing:
            self.lattice = self.createLattice(spacing, boxDimensions)
        else:
            print("The lattice spacing cannot be less than the particle radius!")
            exit(-1)



    def createLattice(self, spacing, boxDimensions):
        particles = []
        [x,y,z] = boxDimensions

        xSites = int((x/2) // spacing)
        validXs = np.unique(np.linspace((x/2)-xSites*spacing, (x/2)+xSites*spacing, 2*xSites+1)%x)

        ySites = int((y/2) // spacing)
        validYs = np.unique(np.linspace((y/2)-ySites*spacing, (y/2)+ySites*spacing, 2*ySites+1)%y)

        zSites = int((z/2) // spacing)
        validZs = np.unique(np.linspace((z/2)-zSites*spacing, (z/2)+zSites*spacing, 2*zSites+1)%z)

        # Start by smushing our valid lattice sites together
        smushed = np.meshgrid(validXs, validYs, validZs)

        # Resshape handles getting all the combinations
        particles = np.array(smushed).T.reshape(-1, 3)

        return particles

