import numpy as np
import matplotlib.pyplot as plt


class ChollaOnTheFlyAnalysisViz:
    '''
    Cholla OnTheFly Visualization Analysis object
        This object will take Cholla On The Fly Analysis objects and do cool plots
    
        Initialized with:
        - OTFAnalysis (ChollaOnTheFlyAnalysis): Cholla On The Fly Analysis object 
        - ImgDir (str): path of directory to place plots
        - dpi (int): (optional) dots-per-inch to specify plot resolution
        - dtype (datatype): (optional) precision to calculate data 
    '''

    def __init__(self, OTFAnalysis, ImgDir, dpi=256):
        self.OTFAnalysis = OTFAnalysis
        self.ImgDir = ImgDir
        self.dpi = dpi

    def get_phasefPath(self, img_ext='png'):
        '''
        Returns the file path for phase space plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"OTFPhaseSpace_{self.OTFAnalysis.nAnalysis:.0f}" + "." + img_ext


    def get_powerSpectrumfPath(self, img_ext='png'):
        '''
        Returns the file path for Power Spectrum plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''


        return self.ImgDir + "/" + f"OTFPowerSpectrum_{self.OTFAnalysis.nAnalysis:.0f}" + "." + img_ext


    def get_combofPath(self, img_ext='png'):
        '''
        Returns the file path for combined plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''


        return self.ImgDir + "/" + f"OTFCombo_{self.OTFAnalysis.nAnalysis:.0f}" + "." + img_ext


    def plot_phase(self, save_plot=False, img_ext='png'):
        '''
        Plot the phase space for this On-The-Fly Analysis object

        Args:
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''
        
        phaseSpaceObj = self.OTFAnalysis.get_phasespace_obj()
        phasespace = phaseSpaceObj.get_phasespace()

        density_min = phaseSpaceObj.OTFPhaseSpaceHead.density_min
        density_max = phaseSpaceObj.OTFPhaseSpaceHead.density_max

        temp_min = phaseSpaceObj.OTFPhaseSpaceHead.T_min
        temp_max = phaseSpaceObj.OTFPhaseSpaceHead.T_max

        curr_z = self.OTFAnalysis.current_z

        # create fig & ax
        fig, ax = plt.subplots(nrows=1, ncols=1)

        # set extents & vmin/max
        phase_extent = np.array([density_min, density_max, temp_min, temp_max])
        l_phase_extent = np.log10(phase_extent)
        vmin, vmax = -8.5, -3

        # call imshow
        l_phasespace = np.log10(phasespace)
        im0 = ax.imshow(l_phasespace, origin='lower', extent = l_phase_extent,
                        vmin = vmin, vmax = vmax)

        # add x/y labels
        xlabel_str = r"$\log_{10} \Delta$"
        _ = ax.set_xlabel(xlabel_str)

        ylabel_str = r"$\log_{10} \rm{T}\ [K]$"
        _ = ax.set_ylabel(ylabel_str)

        # add  x/y lims
        xlow, xupp = -2, 4
        _ = ax.set_xlim(xlow, xupp)

        ylow, yupp = 2, 8
        _ = ax.set_ylim(ylow, yupp)

        # add phase space colorbar
        cbar_ax = fig.add_axes([0.85, 0.1, 0.015, 0.82])

        _ = fig.colorbar(im0, cax=cbar_ax, orientation="vertical")
        _ = cbar_ax.yaxis.set_ticks_position('right')

        # add colorbar label & ensure no overlap w/ticks
        cbar_str = r"$\log_{10} \rm{P}(\Delta, \rm{T})$"
        _ = cbar_ax.set_ylabel(cbar_str, rotation=270)
        _ = cbar_ax.yaxis.labelpad = 20

        # add title
        title_str = rf"$z = {curr_z:.4f}$"
        _ = ax.set_title(title_str)

        _ = fig.tight_layout()

        if save_plot:
            plt.savefig(self.get_phasefPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()


    def plot_powerspectrum(self, dlogk, save_plot=False, img_ext='png'):
        '''
        Plot the Flux Power Spectrum for this On-The-Fly Analysis object

        Args:
            dlogk (float): differential step in log k-space
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        powerSpectrumObj = self.OTFAnalysis.get_powerspectrum_obj(dlogk)
        kvals = powerSpectrumObj.get_kvals()
        powerspectrum = powerSpectrumObj.get_powerspectrum()
        curr_z = self.OTFAnalysis.current_z

        delta2F = (1./np.pi) * kvals * powerspectrum

        # currently, power spectrum has only 1024 data points that are NOT zero. plot only those
        goodPS_ind = np.argwhere(powerspectrum != 0).flatten()

        kvals_good = kvals[goodPS_ind]
        delta2F_good = delta2F[goodPS_ind]

        fig, ax = plt.subplots(nrows=1, ncols=1)

        _ = ax.plot(kvals_good, delta2F_good, marker='.', markersize=10)

        # add x/y labels
        xlabel_str = r'$k\ [\rm{s\ km^{-1}}] $'
        _ = ax.set_xlabel(xlabel_str)

        ylabel_str = r'$\Delta_F^2 (k)$'
        _ = ax.set_ylabel(ylabel_str)

        # add  x/y lims
        xlow, xupp = 1e-3, 4e-1
        _ = ax.set_xlim(xlow, xupp)

        ylow, yupp = 1e-6, 1e2
        _ = ax.set_ylim(ylow, yupp)

        # set x/y log-scale
        _ = ax.set_xscale('log')
        _ = ax.set_yscale('log')

        # add title
        title_str = rf"$z = {curr_z:.4f}$"
        _ = ax.set_title(title_str)

        # add background grid
        _ = ax.grid(which='both', axis='both', alpha=0.3)

        _ = fig.tight_layout()

        if save_plot:
            plt.savefig(self.get_powerSpectrumfPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()


    def plot_combo(self, dlogk, save_plot=False, img_ext='png'):
        '''
        Plot the Flux Power Spectrum and Phase Space for this On-The-Fly Analysis object

        Args:
            dlogk (float): differential step in log k-space
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        phaseSpaceObj = self.OTFAnalysis.get_phasespace_obj()

        # set phase space info
        phasespace = phaseSpaceObj.get_phasespace()

        density_min = phaseSpaceObj.OTFPhaseSpaceHead.density_min
        density_max = phaseSpaceObj.OTFPhaseSpaceHead.density_max

        temp_min = phaseSpaceObj.OTFPhaseSpaceHead.T_min
        temp_max = phaseSpaceObj.OTFPhaseSpaceHead.T_max

        curr_z = self.OTFAnalysis.current_z

        # set power spectrum info 

        powerSpectrumObj = self.OTFAnalysis.get_powerspectrum_obj(dlogk)
        kvals = powerSpectrumObj.get_kvals()
        powerspectrum = powerSpectrumObj.get_powerspectrum()
        curr_z = self.OTFAnalysis.current_z

        delta2F = (1./np.pi) * kvals * powerspectrum

        goodPS_ind = np.argwhere(powerspectrum != 0).flatten()

        kvals_good = kvals[goodPS_ind]
        delta2F_good = delta2F[goodPS_ind]


        fig, (ax_PS, ax_phase) = plt.subplots(nrows=1, ncols=2, figsize=(10.5,6))

        # set extents & vmin/max
        phase_extent = np.array([density_min, density_max, temp_min, temp_max])
        l_phase_extent = np.log10(phase_extent)
        vmin, vmax = -8.5, -3

        # call imshow and plot
        l_phasespace = np.log10(phasespace)
        im0 = ax_phase.imshow(l_phasespace, origin='lower', extent = l_phase_extent,
                        vmin = vmin, vmax = vmax)
        _ = ax_PS.plot(kvals_good, delta2F_good, marker='.', markersize=10)

        # add x labels
        xlabel_str = r"$\log_{10} \Delta$"
        _ = ax_phase.set_xlabel(xlabel_str)

        xlabel_str = r'$k\ [\rm{s\ km^{-1}}] $'
        _ = ax_PS.set_xlabel(xlabel_str)

        # add y labels
        ylabel_str = r"$\log_{10} \rm{T}\ [K]$"
        _ = ax_phase.set_ylabel(ylabel_str)

        ylabel_str = r'$\Delta_F^2 (k)$'
        _ = ax_PS.set_ylabel(ylabel_str)

        # add x lims
        xlow, xupp = -2, 4
        _ = ax_phase.set_xlim(xlow, xupp)
        
        xlow, xupp = 1e-3, 4e-1
        _ = ax_PS.set_xlim(xlow, xupp)

        # add y lims
        ylow, yupp = 2, 8
        _ = ax_phase.set_ylim(ylow, yupp)

        ylow, yupp = 1e-6, 1e2
        _ = ax_PS.set_ylim(ylow, yupp)

        # set x/y powerspectrum to log-scale
        _ = ax_PS.set_xscale('log')
        _ = ax_PS.set_yscale('log')

        # add background grid
        _ = ax_PS.grid(which='both', axis='both', alpha=0.3)

        # add phase space colorbar
        _ = fig.subplots_adjust(left=0.1, top=0.9)
        cbar_ax = fig.add_axes([0.98, 0.145, 0.015, 0.699])
        _ = fig.colorbar(im0, cax=cbar_ax, orientation="vertical")
        _ = cbar_ax.yaxis.set_ticks_position('right')

        # add colorbar label & ensure no overlap w/ticks
        cbar_str = r"$\log_{10} \rm{P}(\Delta, \rm{T})$"
        _ = cbar_ax.set_ylabel(cbar_str, rotation=270)
        _ = cbar_ax.yaxis.labelpad = 20

        # add title
        title_str = rf"$z = {curr_z:.4f}$"
        _ = fig.suptitle(title_str, y=0.95)

        _ = fig.tight_layout()

        if save_plot:
            plt.savefig(self.get_combofPath(img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()
