import numpy as np

class Polymer:

    def __init__(self, length, radius):
        self.monomers = [] # Will be filled using initializePolymer() in sim!
        self.numMonomers = length
        self.monomerRadius = radius