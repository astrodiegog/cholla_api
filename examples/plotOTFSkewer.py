import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from cholla_api.run.ChollaGrid import ChollaGrid
from cholla_api.OTFanalysis.ChollaOnTheFlySkewers import ChollaOnTheFlySkewers
from cholla_api.viz.ChollaOnTheFlySkewerViz import ChollaOnTheFlySkewerViz

plt.style.use("dstyle")
plt.figure()

# grab necessary directory paths
baseCosmoPath = str(Path().absolute())
OTFSkewersDir0 = baseCosmoPath + "/skewers_restart0"
OTFSkewersDir1 = baseCosmoPath + "/skewers_restart1"
imgDir0 = baseCosmoPath + "/skewer_plots0"
imgDir1 = baseCosmoPath + "/skewer_plots1"

# create ChollaGrid object
nprocs = 8
nx, ny, nz = 512, 512, 512
xmin, xmax = 0., 50000.0
ymin, ymax = 0., 50000.0
zmin, zmax = 0., 50000.0
grid = ChollaGrid(nprocs, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax)

# set nstride
nStride = 4

# set number of outputs before restart
nFilesPreRestart = 8

# set Skewer Outputs to plot
nSkewerOutputs = np.arange(18)

# set Skewer x-IDs to plot
nSkewerIDs = np.arange(10)

print("Starting plot loop for : ", baseCosmoPath.split('/')[-1])
for nOutput in nSkewerOutputs:
    if nOutput >= nFilesPreRestart:
        OTFSkewersDir = OTFSkewersDir1
        imgDir = imgDir1
        nOutput = int(nOutput - nFilesPreRestart + 1)
    else:
        OTFSkewersDir = OTFSkewersDir0
        imgDir = imgDir0

    # create ChollaOnTheFlySkewers object
    OTFSkewers = ChollaOnTheFlySkewers(nOutput, OTFSkewersDir, grid)

    # grab current redshift
    curr_z = OTFSkewers.current_z

    # grab x-ChollaOTFSkewer object
    OTFSkewersx = OTFSkewers.get_skewersx_obj(nStride)

    # grab length over x-cell
    dx = grid.dx

    for nSkewerID in nSkewerIDs:
        print(f"\t\tPlotting for skew id : {nSkewerID:.0f}")

        # grab ChollaOTFSkewer object
        OTFSkewx = OTFSkewersx.get_skewer_obj(nSkewerID)

        # create ChollaOTFSkewer Viz object
        OTFSkewxViz = ChollaOnTheFlySkewerViz(OTFSkewx, imgDir)

        OTFSkewxViz.plot_combo(dx, curr_z, nSkewerOutput=nOutput, save_plot=True)

    print(f"\tDone plotting for output : {nOutput:.0f}")

print("Completed plot loop for : ", baseCosmoPath.split('/')[-1])

