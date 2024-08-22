import numpy as np
import matplotlib.pyplot as plt


class ChollaOnTheFlySkewerViz:
    '''
    Cholla OnTheFly Visualization Skewer object
        This object will take Cholla On The Fly Skewer objects and do cool plots
    
        Initialized with:
        - OTFSkew (ChollaOnTheFlySkewer): Cholla On The Fly Skewer object 
        - ImgDir (str): path of directory to place plots
        - dpi (int): (optional) dots-per-inch to specify plot resolution
        - dtype (datatype): (optional) precision to calculate data 
    '''

    def __init__(self, OTFSkew, ImgDir, dpi=256):
        self.OTFSkew = OTFSkew
        self.ImgDir = ImgDir
        self.dpi = dpi

    def get_combofPath(self, nSkewerOutput, img_ext='png'):
        '''
        Returns the file path for combined plot

        Args:
	    nSkewerOutput (int): number of the skewer output
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"OTFSkew_nOut{nSkewerOutput:.0f}_{self.OTFSkew.OTFSkewerHead.skew_key[-1]}id{self.OTFSkew.OTFSkewerHead.skew_id:.0f}" + "." + img_ext


    def plot_combo(self, di, z, save_plot=False, nSkewerOutput=None, img_ext='png'):
        '''
        Plot all saved quantities in a skewer calculated On-TheFly

        Args:
            di (float): cell size in i-dimension
            z (float): current redshift of the skewer
            nSkewerOutput (int): (optional) number of the skewer output, used for file name
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

	# coordinate positions along line-of-sight
        los_arr = np.arange(self.OTFSkew.OTFSkewerHead.n_i) * di

        # grab datasets
        density = self.OTFSkew.get_density()
        HIdensity = self.OTFSkew.get_HIdensity()
        HeIIdensity = self.OTFSkew.get_HeIIdensity()
        vel = self.OTFSkew.get_losvelocity()
        temp = self.OTFSkew.get_temperature()

        # define limits for each dataset
        density_scale = -1., 4.
        HIdensity_scale = -8., 0.
        HeIIdensity_scale = -6., 3.
        vel_scale = -350., 350.
        temp_scale = 3.5, 8.

        fig, ax_all = plt.subplots(nrows=5, ncols=1, figsize=(10,8))

        # define axis objects for each dataset
        ax_density = ax_all[0]
        ax_HIdensity = ax_all[1]
        ax_HeIIdensity = ax_all[2]
        ax_vel = ax_all[3]
        ax_temp = ax_all[4]

        # scatter plot the data
        _ = ax_density.scatter(los_arr, np.log10(density))
        _ = ax_HIdensity.scatter(los_arr, np.log10(HIdensity))
        _ = ax_HeIIdensity.scatter(los_arr, np.log10(HeIIdensity))
        _ = ax_vel.scatter(los_arr, vel)
        _ = ax_temp.scatter(los_arr, np.log10(temp))

        # apply y lims
        _ = ax_density.set_ylim(density_scale[0], density_scale[1])
        _ = ax_HIdensity.set_ylim(HIdensity_scale[0], HIdensity_scale[1])
        _ = ax_HeIIdensity.set_ylim(HeIIdensity_scale[0], HeIIdensity_scale[1])
        _ = ax_vel.set_ylim(vel_scale[0], vel_scale[1])
        _ = ax_temp.set_ylim(temp_scale[0], temp_scale[1])

        # remove bottom ticks for all but bottom-most plot
        for ax in ax_all[:-1]:
            _ = ax.tick_params(axis='both', which='both', labelbottom=False)
    
        # add background grid
        for ax in ax_all:
            _ = ax.grid(which='both', axis='both', alpha=0.3)

        # define x-label coordinate as 5% of length of skewer
        xlabel_coord = self.OTFSkew.OTFSkewerHead.n_i * di * 0.05

        # add labels with annotations
        densityunits_str = r'$\log_{10}[h^2 \rm{M}_{\odot} \rm{kpc}^{-3}]$'
        density_label = r'$\rho - $' + densityunits_str
        density_coords = xlabel_coord, 2.5
        _ = ax_density.annotate(density_label, xy=density_coords)

        HIdensity_label = r'$\rho_{\rm{HI}} - $' + densityunits_str
        HIdensity_coords = xlabel_coord, -2.5
        _ = ax_HIdensity.annotate(HIdensity_label, xy=HIdensity_coords)

        HeIIdensity_label = r'$\rho_{\rm{HeII}} - $' + densityunits_str
        HeIIdensity_coords = xlabel_coord, 0.5
        _ = ax_HeIIdensity.annotate(HeIIdensity_label, xy=HeIIdensity_coords)

        vel_label = r'$v_{\rm{LOS}} - \rm{km}\ \rm{s}^{-1}$'
        vel_coords = xlabel_coord, 100.
        _ = ax_vel.annotate(vel_label, xy=vel_coords)
	
        temp_label = r'$\rm{T} - \log_{10} [\rm{K}]$'
        temp_coords = xlabel_coord, 7.
        _ = ax_temp.annotate(temp_label, xy=temp_coords)
	
        # add x label
        xlabel_str = r"$h^{-1} \rm{kpc}$"
        ax_temp.set_xlabel(xlabel_str)

        # add skew id
        skew_label = rf'${self.OTFSkew.OTFSkewerHead.skew_key[-1]}$' 
        skew_label += r'$\rm{skewid}$' 
        skew_label += rf'$- {self.OTFSkew.OTFSkewerHead.skew_id}$'
        skew_coords_y = 2.5
        skew_coords_x = self.OTFSkew.OTFSkewerHead.n_i * di
        if self.OTFSkew.OTFSkewerHead.skew_id < 1e1: # one digit
            skew_coords_x *= 0.86
        elif self.OTFSkew.OTFSkewerHead.skew_id < 1e2: # two digit
            skew_coords_x *= 0.84
        elif self.OTFSkew.OTFSkewerHead.skew_id < 1e3: # three digit
            skew_coords_x *= 0.83
        elif self.OTFSkew.OTFSkewerHead.skew_id < 1e4: # four digit
            skew_coords_x *= 0.81
        elif self.OTFSkew.OTFSkewerHead.skew_id < 1e5: # five digit
            skew_coords_x *= 0.79
        skew_coords = skew_coords_x, skew_coords_y
        _ = ax_density.annotate(skew_label, xy=skew_coords)

        # add title
        title_str = rf"$z = {z:.4f}$"
        _ = fig.suptitle(title_str, y=0.94)
        _ = fig.tight_layout()

        if save_plot:
            assert nSkewerOutput is not None
            plt.savefig(self.get_combofPath(nSkewerOutput, img_ext), dpi=self.dpi, bbox_inches = "tight")
        else:
            plt.show()
                           



