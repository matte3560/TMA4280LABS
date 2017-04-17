#!/usr/bin/python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parse args to get files containing results
parser = argparse.ArgumentParser(
        description='Used to plot poisson problem result stored in CSV files.')
parser.add_argument('grid', type=str,
        help='CSV file containing grid data.')
parser.add_argument('u', type=str,
        help='CSV file containing result matrix data.')
parser.add_argument('output', type=str,
        help='Output filename (without extension).')
parser.add_argument('-e', '--eps', action='store_true',
        help='Output to EPS format (for use in LaTeX document).')
options = parser.parse_args()

# Get filenames
filename_grid = options.grid
filename_u = options.u
filename_output = options.output

# Parse data from files
grid = np.genfromtxt(filename_grid, dtype=float, delimiter=";")
u = np.genfromtxt(filename_u, dtype=float, delimiter=";")
n = grid.shape[0]-1

# Create plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(grid[1:n], grid[1:n])
ax.plot_surface(X, Y, u, cmap=plt.cm.plasma)

# Set display region and labels
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('u')

# Output to file
if options.eps:
    fig.savefig( filename_output + '.eps', format='eps', dpi=1000)
else:
    fig.savefig( filename_output + '.png', format='png', dpi=250)
