import numpy as np
from scipy.special import erf
from time import time

from cholla_api.run.ChollaGrid import ChollaGrid
from cholla_api.OTFanalysis.ChollaOnTheFlyAnalysis import ChollaOnTheFlyAnalysis
from cholla_api.analysis.ChollaCosmoCalculator import ChollaCosmologyHead
from cholla_api.analysis.ChollaSkewersCalc import ChollaSkewerCosmoCalculator
from cholla_api.snap.ChollaSnap import ChollaSnapHead
from cholla_api.OTFanalysis.ChollaOnTheFlySkewers import ChollaOnTheFlySkewers

print('hellllllllo')
# save impt paths
baseCosmoPath = "/data/groups/comp-astro/digarza/Lya_Test/512_box_taucalc"
OTFSkewersDir = baseCosmoPath + "/skewers"


# create ChollaGrid object
nprocs = 8
nx, ny, nz = 512, 512, 512
xmin, xmax = 0., 50000.0
ymin, ymax = 0., 50000.0
zmin, zmax = 0., 50000.0
grid = ChollaGrid(nprocs, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax)

# assume cube grid, such that line-of-sight and cell width are the same, just use x info
nlos, dx = nx, grid.dx

# set the stride number
nStride = 4

# set Analysis Outputs to plot
nSkewerOutputs = np.arange(16,26)

isCosmoHeadSet = False

precision = np.float64

for nOutput in nSkewerOutputs:
    print(f"Starting calculations for output {nOutput:.0f}")
    time0 = time()

    # create ChollaOnTheFlySkewers object & grab id=0 along x axis
    OTFSkewers = ChollaOnTheFlySkewers(nOutput, OTFSkewersDir, grid)
    OTFSkewersx = OTFSkewers.get_skewersx_obj(nStride)
    OTFSkewersy = OTFSkewers.get_skewersy_obj(nStride)
    OTFSkewersz = OTFSkewers.get_skewersz_obj(nStride)

    if (isCosmoHeadSet is False):
        # create cosmo header with first output
        OmegaM, OmegaR = OTFSkewers.Omega_M, OTFSkewers.Omega_R
        OmegaK, OmegaL = OTFSkewers.Omega_K, OTFSkewers.Omega_L
        w0, wa = OTFSkewers.w0, OTFSkewers.wa
        H0 = OTFSkewers.H0
        chCosmoHead = ChollaCosmologyHead(OmegaM, OmegaR, OmegaK, OmegaL, w0, wa, H0)
        isCosmoHeadSet = True

    # create skew cosmo calc object
    snapHead = ChollaSnapHead(nOutput + 1)
    snapHead.a = OTFSkewers.current_a
    skewCosmoCalc = ChollaSkewerCosmoCalculator(snapHead, chCosmoHead, nlos, dx, precision)
    
    # create all fluxes for each dimension
    fluxes_x = np.zeros((OTFSkewersx.OTFSkewersiHead.n_skews, OTFSkewersx.OTFSkewersiHead.n_i))
    fluxes_y = np.zeros((OTFSkewersy.OTFSkewersiHead.n_skews, OTFSkewersy.OTFSkewersiHead.n_i))
    fluxes_z = np.zeros((OTFSkewersz.OTFSkewersiHead.n_skews, OTFSkewersz.OTFSkewersiHead.n_i))

    # loop over each OTFSkeweri object
    all_OTFSkewers = [OTFSkewersx, OTFSkewersy, OTFSkewersz]
    for i, OTFSkewer_i in enumerate(all_OTFSkewers):
        if (i==0):
            fluxes = fluxes_x
            print("\tCalculations starting along x-axis")
        elif (i==1):
            fluxes = fluxes_y
            print("\tCalculations starting along y-axis")
        else:
            fluxes = fluxes_z
            print("\tCalculations starting along z-axis")
        time_1 = time()
        for nSkewerID in range(OTFSkewer_i.OTFSkewersiHead.n_skews):
            OTFSkewer = OTFSkewer_i.get_skewer_obj(nSkewerID)
            vel = OTFSkewer.get_losvelocity(precision)
            densityHI = OTFSkewer.get_HIdensity(precision)
            temp = OTFSkewer.get_temperature(precision)
            taus = skewCosmoCalc.optical_depth_Hydrogen(densityHI, vel, temp, use_forloop=True)
            fluxes[nSkewerID] = np.exp(-taus)
        axistime = time() - time_1
        print(f"\t\tCalculations completed for this axis in {axistime:.5e} secs")

    all_fluxes = np.zeros((3, OTFSkewersx.OTFSkewersiHead.n_skews, OTFSkewersx.OTFSkewersiHead.n_i))
    all_fluxes[0] = fluxes_x
    all_fluxes[1] = fluxes_y
    all_fluxes[2] = fluxes_z

    flux_mean_x = np.mean(fluxes_x)
    flux_mean_y = np.mean(fluxes_y)
    flux_mean_z = np.mean(fluxes_z)
    flux_mean_all = np.mean(all_fluxes)

    timef = time()

    print(f"Mean flux along x-axis is {flux_mean_x:.5e}")
    print(f"Mean flux along y-axis is {flux_mean_y:.5e}")
    print(f"Mean flux along z-axis is {flux_mean_z:.5e}")
    print(f"Mean flux over all dimensions is {flux_mean_all:.5e}")
    print(f"Calculations took {timef - time0:.5e} secs")



