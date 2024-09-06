import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from cholla_api.run.ChollaGrid import ChollaGrid
from cholla_api.OTFanalysis.ChollaOnTheFlyAnalysis import ChollaOnTheFlyAnalysis
from cholla_api.viz.ChollaOnTheFlyAnalysisViz import ChollaOnTheFlyAnalysisViz

plt.style.use("dstyle")
plt.figure()

# grab necessary directory paths
baseCosmoPath = str(Path().absolute())
OTFAnalysisDir = baseCosmoPath + "/analysis"
imgDir = baseCosmoPath + "/analysis_plots"

# create ChollaGrid object
nprocs = 512
nx, ny, nz = 2048, 2048, 2048
xmin, xmax = 0., 50000.0
ymin, ymax = 0., 50000.0
zmin, zmax = 0., 50000.0
grid = ChollaGrid(nprocs, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax)

# set dlogk
dlogk = 0.0001

# set Analysis Outputs to plot
nAnalysisOutputs = np.arange(17)

print("Starting plot loop for : ", baseCosmoPath.split('/')[-1])
for nOutput in nAnalysisOutputs:
    # create ChollaOnTheFlyAnalysis object
    OTFAnalysis = ChollaOnTheFlyAnalysis(nOutput, OTFAnalysisDir, grid)
    
    # create ChollaOnTheFlyAnalysisViz object
    viz = ChollaOnTheFlyAnalysisViz(OTFAnalysis, imgDir)

    # call methods to plot phase space, power spectrum, and both together
    viz.plot_phase(save_plot=True)

    viz.plot_powerspectrum(dlogk, save_plot=True)

    viz.plot_combo(dlogk, save_plot=True)

    print(f"\tDone plotting for output : {nOutput:.0f}")

print("Completed plot loop for : ", baseCosmoPath.split('/')[-1])

