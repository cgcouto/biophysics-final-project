# Code heavily based on this tutorial found here: 
# https://aaronolsen.github.io/tutorials/3d_visualization/plot_points.html

# Also the reference manual might be helpful here:
# https://cran.r-project.org/web/packages/svgViewR/svgViewR.pdf

# Load the svgViewR package
library(svgViewR)
library(reticulate)

# Change this to generate a particular animation (npy file minus the file extension)
npyName <- "ecoli_stiffness_emphasized"

# Import numpy and load the given npy filename
np <- import("numpy")
mat <- np$load(paste(npyName, ".npy", sep=''))

# Create a permutation vector
p <- c(2, 3, 1)

# Reshape the matrix using aperm
mat_reshaped <- aperm(mat, p)

# Set number of iterations
n_iter <- 100

# Create animated point array
points3da <- mat_reshaped

# Open a connection to .html file
svg.new(file=paste(npyName, ".html", sep=''))

# Add points and lines to file
svg.points(points3da, col.fill='blue', col.stroke='blue')
svg.lines(points3da)

# Add coordinate axis planes around the points
svg_frame <- svg.frame(points3da)

# Close the file connection
svg.close()

