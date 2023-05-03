import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.spatial.distance import pdist
from util import find_bond_angle, get_positions

# Script generates histograms of bond angles and bond lengths
# for simulation data stored in an npy file 

# Npy file name (including .npy file extension)
npyName = "ecoli_stiffness_emphasized.npy"

frameData = np.load(npyName)

endingMonomers = frameData[-1]

# Get bond angles
groups = np.column_stack((endingMonomers[:-2], endingMonomers[1:-1], endingMonomers[2:]))
angles = np.apply_along_axis(find_bond_angle, 1, groups)

plt.hist(angles, bins=10)
plt.xlim(0, pi)
plt.title("Histogram of bond angles")
plt.xlabel("Bond angle (rad)")
plt.ylabel("Number of occurrences")
plt.xticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi], ["0", "π/6", "π/3", "π/2", "2π/3", "5π/6", "π"])
plt.show()

# Get pairwise distances
pairwiseDists = pdist(endingMonomers)

# Separate bonded vs. non-bonded distances
bondedInds = get_positions(len(endingMonomers))
unbondedInds = np.array(list(set(range(len(pairwiseDists))) - set(bondedInds)))

# Use the inds to grab the relevant bonded dists
bondedDists = pairwiseDists[bondedInds]

plt.hist(bondedDists, bins=10)
plt.xlim(0, 2)
plt.title("Histogram of bond distances")
plt.xlabel("Bond distance")
plt.ylabel("Number of occurrences")
plt.show()