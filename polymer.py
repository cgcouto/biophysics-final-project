import numpy as np

class Polymer:

    def __init__(self, length, diameter):
        self.monomers = [] # Will be filled using initializePolymer() in sim!
        self.numMonomers = length
        self.monomerDiam = diameter