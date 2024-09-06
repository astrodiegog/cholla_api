import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

from cholla_api.run.ChollaGrid import ChollaGrid
from cholla_api.run.ChollaGlobal import ChollaGlobal
from cholla_api.run.ChollaMacroFlags import ChollaMacroFlags
from cholla_api.viz.ChollaCosmoVizAnalysis import ChollaCosmoVizAnalysis
from cholla_api.OTFanalysis.ChollaOnTheFlyAnalysis import ChollaOnTheFlyAnalysis
from cholla_api.analysis.ChollaCosmoGridCalculator import ChollaCosmologyHead

plt.style.use("dstyle")
_ = plt.figure()

# check args for what to plot
args = sys.argv
temp_key, density_key, densityDM_key = 'temp', 'density', 'densityDM'
plottemp_bool = temp_key in args
plotdensity_bool = density_key in args
plotDM_bool = densityDM_key in args


if not (plottemp_bool or plotdensity_bool or plotDM_bool):
    print("oops, no plot specified, exiting!")
    exit()


# grab necessary directory paths
baseCosmoPath = str(Path().absolute())
OTFAnalysisDir = baseCosmoPath + "/analysis"
dataPath = baseCosmoPath + "/data"
imgDir = baseCosmoPath + "/cosmoviz_plots"

# create ChollaGrid object
nprocs = 8
nx, ny, nz = 512, 512, 512
xmin, xmax = 0., 50000.0
ymin, ymax = 0., 50000.0
zmin, zmax = 0., 50000.0
grid = ChollaGrid(nprocs, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax)

# create ChollaMacroFlags object

# set compiler flags
DE_flag = True
Grav_flag, Part_flag, Cosmo_flag = True, True, True
CoolingGrack_flag = False
OTFAnalysis_flag = True
CoolingGPU_flag = True
chMacroFlags = ChollaMacroFlags(DE_flag, Grav_flag, Part_flag, Cosmo_flag, CoolingGrack_flag,
                                    CoolingGPU_flag, OTFAnalysis_flag)

# create ChollaOnTheFlyAnalysis object & grab cosmo params
nOutput = 0
OTFAnalysis = ChollaOnTheFlyAnalysis(nOutput, OTFAnalysisDir, grid)
OmegaM, OmegaR = OTFAnalysis.Omega_M, OTFAnalysis.Omega_R
OmegaK, OmegaL = OTFAnalysis.Omega_K, OTFAnalysis.Omega_L
w0, wa = OTFAnalysis.w0, OTFAnalysis.wa
H0 = OTFAnalysis.H0
chCosmoHead = ChollaCosmologyHead(OmegaM, OmegaR, OmegaK, OmegaL, w0, wa, H0)

# create ChollaGlobal object
chGlobal = ChollaGlobal(baseCosmoPath, grid, chMacroFlags)

# set global attributes with snap=1 box=0
chGlobal.set_globalattrs(1,0,all_units=True)

# grab impt units & set mu
velocity_unit = 1.e5 # [km s-1] --> [cm s-1]
energy_unit = velocity_unit * velocity_unit
gamma = chGlobal.gamma
mu = 0.60

nSnaps = np.arange(1,20)
boxAttrs = 0
for nSnap in nSnaps:
    print(f"--- Starting cosmoviz for snapshot : {nSnap:.0f} ---")
    # create ChollaSnap object
    snap = chGlobal.get_snap(nSnap, box_for_attrs=boxAttrs)

    # create cosmoviz object
    chCosmoViz = ChollaCosmoVizAnalysis(snap, grid, chMacroFlags, chCosmoHead, imgDir)

    if plotdensity_bool:
        # plot gas density diagnostic
        print("density!")
        chCosmoViz.plot_density(save_plot=True)

    if plotDM_bool:
        # plot dark matter density diagnostic
        print("DM!")
        chCosmoViz.plot_densityDM(save_plot=True)

    if plottemp_bool:
        # plot temperature diagnostic
        print("temp!")
        chCosmoViz.plot_temp(gamma, mu, energy_unit, save_plot=True)
