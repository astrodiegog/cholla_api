import numpy as np
import matplotlib.pyplot as plt

from cholla_api.data.ChollaBox import ChollaBox
from cholla_api.data.ChollaBoxHydroCalc import ChollaBoxHydroCalc
from cholla_api.analysis.ChollaCosmoGridCalculator import ChollaCosmoGridCalculator

class ChollaCosmoVizAnalysis:
    '''
    Cholla Cosmological Visualization Analysis object
        This object will create snapshot-wide visualization plots

        Initialized with:
        - ChollaSnap (ChollaSnap): ChollaSnap object
        - ChollaGrid (ChollaGrid): ChollaGrid object
        - ImgDir (str): path of directory to place plots
        - OmegaM (float): present-day energy density parameter for matter
        - OmegaR (float): present-day energy density parameter for radiation
        - OmegaK (float): present-day energy density parameter for spatial curvature
        - OmegaL (float): present-day energy density parameter for dark energy
        - w0 (float): constant term in dark energy equation of state
        - wa (float): linear term in dark energy equation of state
        - H0 (float): present-day Hubble parameter in units of [km / s / Mpc]
        - dpi (int): (optional) dots-per-inch to specify plot resolution
        - dtype (datatype): (optional) precision to calculate data 
    '''
    
    def __init__(self, ChollaSnap, ChollaGrid, ImgDir, OmegaM, OmegaR, OmegaK, OmegaL, w0, wa, H0, dpi=256, precision=np.float32):
        self.snap = ChollaSnap
        self.ImgDir = ImgDir
        
        self.cosmoCalc = ChollaCosmoGridCalculator(ChollaGrid, OmegaM, OmegaR, OmegaK, OmegaL, w0, wa, H0)
        self.boxHeads = ChollaGrid.get_BoxHeads()
        self.nx_global = ChollaGrid.nx_global
        self.ny_global = ChollaGrid.ny_global
        self.nz_global = ChollaGrid.nz_global

        self.xmin = ChollaGrid.xmin
        self.xmax = ChollaGrid.Lx + self.xmin

        self.ymin = ChollaGrid.ymin
        self.ymax = ChollaGrid.Ly + self.ymin

        self.zmin = ChollaGrid.zmin
        self.zmax = ChollaGrid.Lz + self.zmin

        self.precision = precision
        self.dpi = dpi

        # number of histogram bins
        self.nhist_bins = 100


    def get_densityfPath(self, img_ext='png'):
        '''
        Returns the file path for density plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"densitydiagnostic_{self.snap.SnapHead.nSnap:.0f}" + "." + img_ext

    
    def get_densityDMfPath(self, img_ext='png'):
        '''
        Returns the file path for dark matter density plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"densityDMdiagnostic_{self.snap.SnapHead.nSnap:.0f}" + "." + img_ext


    def get_tempfPath(self, img_ext='png'):
        '''
        Returns the file path for temperature plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"tempdiagnostic_{self.snap.SnapHead.nSnap:.0f}" + "." + img_ext


    def plot_density(self, save_plot=False, img_ext='png'):
        '''
        Plot the density diagnostic for this snapshot.
            Create a 2x2 plot:
                top-left will be histogram distribution of all density
                top-right will be density projection along x-axis
                bottom-left will be density projection along y-axis
                bottom-right will be density projection along z-axis
        
        Orientation is such that the "inner" axes will be shared among the three
            projection plots. Two bottom plots share x-axis. Two right plots
            share y-axis.

        imshow() uses first axis as the vertical axis
        n_x (y,z) needs to be transposed so vertical is z and horizontal is y
        n_z (x,y) needs to flip order of vertical so x is decreasing downward
        n_y (x,z) needs to flip order of both axes so x is decreasing downward 
            and z is decreasing leftward

        Assumes density is saved in [h2 Msun / kpc3] comoving density units (!)

        Args:
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        print("Starting density plot calcs")

        # in units of [h2 Msun / kpc3]
        l_histmin, l_histmax = -3.0, 3.0
        l_cbarmin, l_cbarmax = 1.0, 5.0
        # set the top-left histogram bins
        l_nhist_bins = np.linspace(l_histmin, l_histmax, self.nhist_bins)
        
        # set four arrays to be plotted in log-space
        l_nhist_snap = np.zeros(self.nhist_bins - 1, dtype=self.precision)
        nx_snap = np.zeros((self.ny_global, self.nz_global), dtype=self.precision)
        ny_snap = np.zeros((self.nx_global, self.nz_global), dtype=self.precision)
        nz_snap = np.zeros((self.nx_global, self.ny_global), dtype=self.precision)

        for boxhead in self.boxHeads:
            print(f"\tStarting calcs for boxhead : {boxhead.nBox:.0f}")
            box = ChollaBox(self.snap.SnapPath, boxhead)
            boxhydrocalc = ChollaBoxHydroCalc(box, self.precision)

            # grab density
            density = box.get_hydrodata(box.density_str, dtype=self.precision)

            # histogram
            l_nhist_box, _ = np.histogram(np.log10(density.flatten()), bins=l_nhist_bins)

            # calculate xyz projections
            nx_box = boxhydrocalc.Calculator.densityk_projection(density, 0)
            ny_box = boxhydrocalc.Calculator.densityk_projection(density, 1)
            nz_box = boxhydrocalc.Calculator.densityk_projection(density, 2)

            # grab offsets to place onto global array
            startX, startY, startZ = box.BoxHead.offset
            localX, localY, localZ = box.BoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

            # place onto global array
            l_nhist_snap += l_nhist_box
            nx_snap[startY:endY, startZ:endZ] += nx_box
            ny_snap[startX:endX, startZ:endZ] += ny_box
            nz_snap[startX:endX, startY:endY] += nz_box
            print(f"\t\tCompleted calcs for boxhead : {boxhead.nBox:.0f}")
   
        # take log of projections
        l_nx_snap = np.log10(nx_snap)
        l_ny_snap = np.log10(ny_snap)
        l_nz_snap = np.log10(nz_snap)

        # normalize histogram to become a PDF
        totnumcells = self.nx_global * self.ny_global * self.nz_global
        l_nhist_snap_norm = l_nhist_snap / (totnumcells)

        print("Starting plotting !")

        # time to plot (!)
        fig, ax = plt.subplots(figsize=(10,10), nrows=2, ncols=2)
        
        # call imshow and bar
        centers = (l_nhist_bins[:-1] + l_nhist_bins[1:]) / 2.
        width = np.diff(l_nhist_bins)

        im0 = ax[0,0].bar(centers, l_nhist_snap_norm, width=width, align='center', facecolor='C0', edgecolor='none')

        im1 = ax[0,1].imshow(l_nx_snap.T, origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.zmin, self.zmax])

        im2 = ax[1,0].imshow(l_ny_snap[::-1,::-1], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.zmax, self.zmin, self.xmax, self.xmin])

        im3 = ax[1,1].imshow(l_nz_snap[::-1,:], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.xmax, self.xmin])


        # get rid of labels on projection plots
        _ = ax[1,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[0,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[1,0].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)


        # place ticks on left and top
        _ = ax[1,1].tick_params(axis='both', which='both',
                                left=True, top=True, bottom=False, right=False,
                                color='white')

        # place ticks on left and bottom
        _ = ax[0,1].tick_params(axis='both', which='both',
                                left=True, top=False, bottom=True, right=False,
                                color='white')

        # place ticks on right and top
        _ = ax[1,0].tick_params(axis='both', which='both',
                                left=False, top=True, bottom=False, right=True,
                                color='white')

        # add background grid
        _ = ax[0,0].grid(which='both', axis='both', alpha=0.3)

        # set annotating text to label each projection
        rhox_str = r'$\rho_x$'
        rhoy_str = r'$\rho_y$'
        rhoz_str = r'$\rho_z$'
        rhox_str_coords = self.ymin + 0.05*self.ymax, self.zmin + 0.9*self.zmax
        rhoy_str_coords = self.zmin + 0.95*self.zmax, self.xmin + 0.1*self.xmax
        rhoz_str_coords = self.ymin + 0.05*self.ymax, self.xmin + 0.1*self.xmax
        _ = ax[0,1].annotate(rhox_str, xy=rhox_str_coords, color='white', fontsize=35)
        _ = ax[1,0].annotate(rhoy_str, xy=rhoy_str_coords, color='white', fontsize=35)
        _ = ax[1,1].annotate(rhoz_str, xy=rhoz_str_coords, color='white', fontsize=35)

        # set scale
        perc_Lbox = 0.2 # percentage of total box length to use
        scaleunits_str = r"$h^{-1} \rm{Mpc}$"
        scalex_str = rf"${perc_Lbox * (self.xmax - self.xmin) / 1e3 :.0f} \ $" + scaleunits_str
        scaley_str = rf"${perc_Lbox * (self.ymax - self.ymin) / 1e3 :.0f} \ $" + scaleunits_str
        scalez_str = rf"${perc_Lbox * (self.zmax - self.zmin) / 1e3 :.0f} \ $" + scaleunits_str
        scalex_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.zmin + 0.12*self.zmax
        scaley_str_coords = (self.zmax - self.zmin) * (0.95 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax
        scalez_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax

        # place length scale + labels
        lw_Lbox = 3 # linewidth to use
        _ = ax[0,1].axhline(y = self.zmin + 0.1*self.zmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,0].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,1].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[0,1].annotate(scalex_str, xy=scalex_str_coords, color='white', fontsize=20)
        _ = ax[1,0].annotate(scaley_str, xy=scaley_str_coords, color='white', fontsize=20)
        _ = ax[1,1].annotate(scalez_str, xy=scalez_str_coords, color='white', fontsize=20)

        # add shared colorbar to right
        fig.subplots_adjust(right=0.9, bottom=0.05)
        cbar_ax = fig.add_axes([1., 0.035, 0.05, 0.89])
        _ = fig.colorbar(im3, cax=cbar_ax)
        _ = cbar_ax.yaxis.set_ticks_position('right')

        # add colorbar label & ensure no overlap w/ticks
        cbar_str = r"$\log_{10} \rm{Density} \ $" + r"$ [h^2 \rm{M}_{\odot} \rm{kpc}^{-3}]$"
        _ = cbar_ax.set_ylabel(cbar_str, rotation=270, fontsize=25)
        _ = cbar_ax.yaxis.labelpad = 30

        # add histogram unit label 
        # _ = ax[0,0].set_xlabel(cbar_str, fontsize=25)

        # add title
        title_str = r"$\rm{Gas \ Density - }$" + rf"$z = {self.snap.SnapHead.z:.4f}$"
        _ = fig.suptitle(title_str)

        _ = fig.tight_layout()

        print("Done plotting!")

        if save_plot:
            plt.savefig(self.get_densityfPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()


    def plot_temp(self, gamma, mu, energy_units, save_plot=False, img_ext='png'):
        '''
        Plot the temperature diagnostic for this snapshot. Assumes constant mean
            molecular weight and ratio of specific heats in all cells.
            Create a 2x2 plot:
                top-left will be histogram distribution of all temperatures
                top-right will be temperature projection along x-axis
                bottom-left will be temperature projection along y-axis
                bottom-right will be temperature projection along z-axis
        
        Orientation is such that the "inner" axes will be shared among the three
            projection plots. Two bottom plots share x-axis. Two right plots
            share y-axis.

        imshow() uses first axis as the vertical axis
        T_x (y,z) needs to be transposed so vertical is z and horizontal is y
        T_z (x,y) needs to flip order of vertical so x is decreasing downward
        T_y (x,z) needs to flip order of both axes so x is decreasing downward 
            and z is decreasing leftward

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            energy_unit (float): unit conversion from energy code units to cgs
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        code2cgs = energy_units

        print("Starting temp plot calcs")

        # in units of [Kelvin]
        l_histmin, l_histmax = 3.0, 7.0
        l_cbarmin, l_cbarmax = 7.0, 10.0
        # set the top-left histogram bins
        l_Thist_bins = np.linspace(l_histmin, l_histmax, self.nhist_bins)

        # set four arrays to be plotted in log-space
        l_Thist_snap = np.zeros(self.nhist_bins - 1, dtype=self.precision)
        Tx_snap = np.zeros((self.ny_global, self.nz_global), dtype=self.precision)
        Ty_snap = np.zeros((self.nx_global, self.nz_global), dtype=self.precision)
        Tz_snap = np.zeros((self.nx_global, self.ny_global), dtype=self.precision)

        for boxhead in self.boxHeads:
            print(f"\tStarting calcs for boxhead : {boxhead.nBox:.0f}")
            box = ChollaBox(self.snap.SnapPath, boxhead)
            boxhydrocalc = ChollaBoxHydroCalc(box, self.precision)

            # grab density and gas energy, keep in code units
            density = box.get_hydrodata(box.density_str, dtype=self.precision)
            gas_energy = box.get_hydrodata(box.gasenergy_str, dtype=self.precision)

            # calculate temperature
            temp = boxhydrocalc.Calculator.gas_temp(gas_energy, density, gamma, 
                                                    mu, code2cgs)

            # histogram
            l_Thist_box, _ = np.histogram(np.log10(temp.flatten()), bins=l_Thist_bins)

            # calculate xyz projections, can reuse densityk_projection
            Tx_box = boxhydrocalc.Calculator.densityk_projection(temp, 0)
            Ty_box = boxhydrocalc.Calculator.densityk_projection(temp, 1)
            Tz_box = boxhydrocalc.Calculator.densityk_projection(temp, 2)

            # grab offsets to place onto global array
            startX, startY, startZ = box.BoxHead.offset
            localX, localY, localZ = box.BoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

            # place onto global array
            l_Thist_snap += l_Thist_box
            Tx_snap[startY:endY, startZ:endZ] += Tx_box
            Ty_snap[startX:endX, startZ:endZ] += Ty_box
            Tz_snap[startX:endX, startY:endY] += Tz_box
            print(f"\t\tCompleted calcs for boxhead : {boxhead.nBox:.0f}")


        # take log of projections
        l_Tx_snap = np.log10(Tx_snap)
        l_Ty_snap = np.log10(Ty_snap)
        l_Tz_snap = np.log10(Tz_snap)

        # normalize histogram to become a PDF
        totnumcells = self.nx_global * self.ny_global * self.nz_global
        l_Thist_snap_norm = l_Thist_snap / (totnumcells)

        print("Starting plotting !")

        # time to plot (!)
        fig, ax = plt.subplots(figsize=(10,10), nrows=2, ncols=2)

        # call imshow and bar
        centers = (l_Thist_bins[:-1] + l_Thist_bins[1:]) / 2.
        width = np.diff(l_Thist_bins)

        im0 = ax[0,0].bar(centers, l_Thist_snap_norm, width=width, align='center', facecolor='C0', edgecolor='none')

        im1 = ax[0,1].imshow(l_Tx_snap.T, origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.zmin, self.zmax])

        im2 = ax[1,0].imshow(l_Ty_snap[::-1,::-1], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.zmax, self.zmin, self.xmax, self.xmin])

        im3 = ax[1,1].imshow(l_Tz_snap[::-1,:], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.xmax, self.xmin])

        
        # get rid of labels on projection plots
        _ = ax[1,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[0,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[1,0].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)


        # place ticks on left and top
        _ = ax[1,1].tick_params(axis='both', which='both',
                                left=True, top=True, bottom=False, right=False,
                                color='white')

        # place ticks on left and bottom
        _ = ax[0,1].tick_params(axis='both', which='both',
                                left=True, top=False, bottom=True, right=False,
                                color='white')

        # place ticks on right and top
        _ = ax[1,0].tick_params(axis='both', which='both',
                                left=False, top=True, bottom=False, right=True,
                                color='white')

        # add background grid
        _ = ax[0,0].grid(which='both', axis='both', alpha=0.3)

        # set annotating text to label each projection
        Tx_str = r'$T_x$'
        Ty_str = r'$T_y$'
        Tz_str = r'$T_z$'
        Tx_str_coords = self.ymin + 0.05*self.ymax, self.zmin + 0.9*self.zmax
        Ty_str_coords = self.zmin + 0.95*self.zmax, self.xmin + 0.1*self.xmax
        Tz_str_coords = self.ymin + 0.05*self.ymax, self.xmin + 0.1*self.xmax
        _ = ax[0,1].annotate(Tx_str, xy=Tx_str_coords, color='white', fontsize=35)
        _ = ax[1,0].annotate(Ty_str, xy=Ty_str_coords, color='white', fontsize=35)
        _ = ax[1,1].annotate(Tz_str, xy=Tz_str_coords, color='white', fontsize=35)

        # set scale
        perc_Lbox = 0.2 # percentage of total box length to use
        scaleunits_str = r"$h^{-1} \rm{Mpc}$"
        scalex_str = rf"${perc_Lbox * (self.xmax - self.xmin) / 1e3 :.0f} \ $" + scaleunits_str
        scaley_str = rf"${perc_Lbox * (self.ymax - self.ymin) / 1e3 :.0f} \ $" + scaleunits_str
        scalez_str = rf"${perc_Lbox * (self.zmax - self.zmin) / 1e3 :.0f} \ $" + scaleunits_str
        scalex_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.zmin + 0.12*self.zmax
        scaley_str_coords = (self.zmax - self.zmin) * (0.95 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax
        scalez_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax

        # place length scale + labels
        lw_Lbox = 3 # linewidth to use
        _ = ax[0,1].axhline(y = self.zmin + 0.1*self.zmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,0].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,1].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[0,1].annotate(scalex_str, xy=scalex_str_coords, color='white', fontsize=20)
        _ = ax[1,0].annotate(scaley_str, xy=scaley_str_coords, color='white', fontsize=20)
        _ = ax[1,1].annotate(scalez_str, xy=scalez_str_coords, color='white', fontsize=20)

        # add shared colorbar to right
        fig.subplots_adjust(right=0.9, bottom=0.05)
        cbar_ax = fig.add_axes([1., 0.035, 0.05, 0.89])
        _ = fig.colorbar(im3, cax=cbar_ax)
        _ = cbar_ax.yaxis.set_ticks_position('right')

        # add colorbar label & ensure no overlap w/ticks
        cbar_str = r"$\log_{10} \rm{Temperature} \ $" + r"$ [\rm{K}]$"
        _ = cbar_ax.set_ylabel(cbar_str, rotation=270, fontsize=25)
        _ = cbar_ax.yaxis.labelpad = 30

        # add title
        title_str = r"$\rm{Gas \ Temperature - }$" + rf"$z = {self.snap.SnapHead.z:.4f}$"
        _ = fig.suptitle(title_str)

        _ = fig.tight_layout()

        print("Done plotting!")

        if save_plot:
            plt.savefig(self.get_tempfPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    
    def plot_densityDM(self, save_plot=False, img_ext='png'):
        '''
        Plot the dark matter density diagnostic for this snapshot using the
            Cloud-in-Cell density saved in particle data files
            Create a 2x2 plot:
                top-left will be histogram distribution of all density
                top-right will be density projection along x-axis
                bottom-left will be density projection along y-axis
                bottom-right will be density projection along z-axis

        Orientation is such that the "inner" axes will be shared among the three
            projection plots. Two bottom plots share x-axis. Two right plots
            share y-axis.

        Assumes density is saved in [h2 Msun / kpc3] comoving density units (!)

        imshow() uses first axis as the vertical axis
        n_x (y,z) needs to be transposed so vertical is z and horizontal is y
        n_z (x,y) needs to flip order of vertical so x is decreasing downward
        n_y (x,z) needs to flip order of both axes so x is decreasing downward
            and z is decreasing leftward

        Args:
            density_unit (float): unit conversion from density code units to cgs
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        print("Starting dark matter density plot calcs")

        # in units of [h2 Msun / kpc3]
        l_histmin, l_histmax = -3.0, 3.0
        l_cbarmin, l_cbarmax = 2.0, 6.0
        # set the top-left histogram bins
        l_nhist_bins = np.linspace(l_histmin, l_histmax, self.nhist_bins)

        # set four arrays to be plotted in log-space
        l_nhist_snap = np.zeros(self.nhist_bins - 1, dtype=self.precision)
        nx_snap = np.zeros((self.ny_global, self.nz_global), dtype=self.precision)
        ny_snap = np.zeros((self.nx_global, self.nz_global), dtype=self.precision)
        nz_snap = np.zeros((self.nx_global, self.ny_global), dtype=self.precision)

        for boxhead in self.boxHeads:
            print(f"\tStarting calcs for boxhead : {boxhead.nBox:.0f}")
            box = ChollaBox(self.snap.SnapPath, boxhead)
            boxhydrocalc = ChollaBoxHydroCalc(box, self.precision)

            # grab density
            density = box.get_particledata(box.densityCIC_str, dtype=self.precision)

            # histogram
            l_nhist_box, _ = np.histogram(np.log10(density.flatten()), bins=l_nhist_bins)

            # calculate xyz projections
            nx_box = boxhydrocalc.Calculator.densityk_projection(density, 0)
            ny_box = boxhydrocalc.Calculator.densityk_projection(density, 1)
            nz_box = boxhydrocalc.Calculator.densityk_projection(density, 2)

            # grab offsets to place onto global array
            startX, startY, startZ = box.BoxHead.offset
            localX, localY, localZ = box.BoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

            # place onto global array
            l_nhist_snap += l_nhist_box
            nx_snap[startY:endY, startZ:endZ] += nx_box
            ny_snap[startX:endX, startZ:endZ] += ny_box
            nz_snap[startX:endX, startY:endY] += nz_box
            print(f"\t\tCompleted calcs for boxhead : {boxhead.nBox:.0f}")
        

        # take log of projections
        l_nx_snap = np.log10(nx_snap)
        l_ny_snap = np.log10(ny_snap)
        l_nz_snap = np.log10(nz_snap)

        # normalize histogram to become a PDF
        totnumcells = self.nx_global * self.ny_global * self.nz_global
        l_nhist_snap_norm = l_nhist_snap / (totnumcells)

        print("Starting plotting !")

        # time to plot (!)
        fig, ax = plt.subplots(figsize=(10,10), nrows=2, ncols=2)

        # call imshow and bar
        centers = (l_nhist_bins[:-1] + l_nhist_bins[1:]) / 2.
        width = np.diff(l_nhist_bins)

        im0 = ax[0,0].bar(centers, l_nhist_snap_norm, width=width, align='center', facecolor='C0', edgecolor='none')

        im1 = ax[0,1].imshow(l_nx_snap.T, origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.zmin, self.zmax])

        im2 = ax[1,0].imshow(l_ny_snap[::-1,::-1], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.zmax, self.zmin, self.xmax, self.xmin])

        im3 = ax[1,1].imshow(l_nz_snap[::-1,:], origin='lower', vmin=l_cbarmin, vmax=l_cbarmax,
                             extent=[self.ymin, self.ymax, self.xmax, self.xmin])


        # get rid of labels on projection plots
        _ = ax[1,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[0,1].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)
        _ = ax[1,0].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)


        # place ticks on left and top
        _ = ax[1,1].tick_params(axis='both', which='both',
                                left=True, top=True, bottom=False, right=False,
                                color='white')

        # place ticks on left and bottom
        _ = ax[0,1].tick_params(axis='both', which='both',
                                left=True, top=False, bottom=True, right=False,
                                color='white')

        # place ticks on right and top
        _ = ax[1,0].tick_params(axis='both', which='both',
                                left=False, top=True, bottom=False, right=True,
                                color='white')

        # add background grid
        _ = ax[0,0].grid(which='both', axis='both', alpha=0.3)

        # set annotating text to label each projection
        rhox_str = r'$\rho_x$'
        rhoy_str = r'$\rho_y$'
        rhoz_str = r'$\rho_z$'
        rhox_str_coords = self.ymin + 0.05*self.ymax, self.zmin + 0.9*self.zmax
        rhoy_str_coords = self.zmin + 0.95*self.zmax, self.xmin + 0.1*self.xmax
        rhoz_str_coords = self.ymin + 0.05*self.ymax, self.xmin + 0.1*self.xmax
        _ = ax[0,1].annotate(rhox_str, xy=rhox_str_coords, color='white', fontsize=35)
        _ = ax[1,0].annotate(rhoy_str, xy=rhoy_str_coords, color='white', fontsize=35)
        _ = ax[1,1].annotate(rhoz_str, xy=rhoz_str_coords, color='white', fontsize=35)

        # set scale
        perc_Lbox = 0.2 # percentage of total box length to use
        scaleunits_str = r"$h^{-1} \rm{Mpc}$"
        scalex_str = rf"${perc_Lbox * (self.xmax - self.xmin) / 1e3 :.0f} \ $" + scaleunits_str
        scaley_str = rf"${perc_Lbox * (self.ymax - self.ymin) / 1e3 :.0f} \ $" + scaleunits_str
        scalez_str = rf"${perc_Lbox * (self.zmax - self.zmin) / 1e3 :.0f} \ $" + scaleunits_str
        scalex_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.zmin + 0.12*self.zmax
        scaley_str_coords = (self.zmax - self.zmin) * (0.95 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax
        scalez_str_coords = (self.ymax - self.ymin) * (0.05 + 0.*perc_Lbox), self.xmin + 0.88*self.xmax

        # place length scale + labels
        lw_Lbox = 3 # linewidth to use
        _ = ax[0,1].axhline(y = self.zmin + 0.1*self.zmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,0].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[1,1].axhline(y = self.xmin + 0.9*self.xmax, xmin=0.1, xmax=0.1 + perc_Lbox, color='white', lw=lw_Lbox)
        _ = ax[0,1].annotate(scalex_str, xy=scalex_str_coords, color='white', fontsize=20)
        _ = ax[1,0].annotate(scaley_str, xy=scaley_str_coords, color='white', fontsize=20)
        _ = ax[1,1].annotate(scalez_str, xy=scalez_str_coords, color='white', fontsize=20)

        # add shared colorbar to right
        fig.subplots_adjust(right=0.9, bottom=0.05)
        cbar_ax = fig.add_axes([1., 0.035, 0.05, 0.89])
        _ = fig.colorbar(im3, cax=cbar_ax)
        _ = cbar_ax.yaxis.set_ticks_position('right')

        # add colorbar label & ensure no overlap w/ticks
        cbar_str = r"$\log_{10} \rm{Density} \ $" + r"$ [h^2 \rm{M}_{\odot} \rm{kpc}^{-3}]$"
        _ = cbar_ax.set_ylabel(cbar_str, rotation=270, fontsize=25)
        _ = cbar_ax.yaxis.labelpad = 30

        # add title
        title_str = r"$\rm{Dark \ Matter \ Density - }$" + rf"$z = {self.snap.SnapHead.z:.4f}$"
        _ = fig.suptitle(title_str)

        _ = fig.tight_layout()

        print("Done plotting!")

        if save_plot:
            plt.savefig(self.get_densityDMfPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()


