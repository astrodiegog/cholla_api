import numpy as np
import matplotlib.pyplot as plt


class ChollaVizAnalysis:
    '''
    Cholla Visualization Analysis object
        This object will take any data and place it onto a cool output plot
    '''
    
    def __init__(self, dims, offset, dpi=256):
        self.xmin, self.ymin, self.zmin = offset
        self.xdim, self.ydim, self.zdim = dims
        self.xmax = self.xmin + self.xdim
        self.ymax = self.ymin + self.ydim
        self.zmax = self.zmin + self.zdim
        self.dpi = dpi
        
    
    def cosmo_diagnostic(self, phasespace, phase_xedge, phase_yedge, n_x, n_y, n_z, show_ticks=True, fname=None, time=None):
        '''
        Create a 2x2 plot of n_i and phase space
            top-left will be phase space diagram
            top-right will be projection along x-axis
            bottom-left will be projection along y-axis
            bottom-right will be projection along z-axis
        
        Orientation is such that the "inner" axes will be shared among the three 
            projection plots. Two bottom plots share x-axis. Two right plots 
            share y-axis.

        imshow() uses first axis as the vertical axis
        n_x (y,z) needs to be transposed so vertical is z and horizontal is y
        n_z (x,y) needs to flip order of vertical so x is decreasing downward
        n_y (x,z) needs to flip order of both axes so x is decreasing downward 
            and z is decreasing leftward
        
        time is tuple (snap number, cosmic age, redshift)
        '''
        
        # normalize 2D histogram + take its log10
        hist_norm = phasespace / np.sum(phasespace)
        log_histnorm = np.log10(hist_norm)
        
        
        fig, ax = plt.subplots(figsize=(10,10), nrows=2, ncols=2)

        # hard set norm to be 1e-8 to 1e-2 in default binning from ChollaAnalysis
        vmin, vmax = -8, -2
        im0 = ax[0,0].imshow(log_histnorm.T, origin='lower', vmin=vmin, vmax=vmax,
                             extent=[phase_xedge[0], phase_xedge[-1], 
                                     phase_yedge[0], phase_yedge[-1]])

        im1 = ax[0,1].imshow(np.log10(n_x.T), origin='lower', 
                             extent=[self.ymin, self.ymax, self.zmin, self.zmax])

        im2 = ax[1,0].imshow(np.log10(n_y[::-1,::-1]), origin='lower', 
                             extent=[self.zmax, self.zmin, self.xmax, self.xmin])

        im3 = ax[1,1].imshow(np.log10(n_z[::-1,:]), origin='lower',
                             extent=[self.ymin, self.ymax, self.xmax, self.xmin])

        # set phase space labels, place density label+ticks on top
        phase_xlbl = r"$\log_{10} \rho / \bar{\rho}$"
        phase_ylbl = r"$\log_{10}$ T/K"
        ax[0,0].set_xlabel(phase_xlbl), ax[0,0].set_ylabel(phase_ylbl)
        ax[0,0].tick_params(axis='x', which='both', 
                            labelbottom=False, labeltop=True, # move label to top
                            left=False, top=True, bottom=False, right=False,  # place ticks on top
                           )
        ax[0,0].xaxis.set_label_position("top")

        # get rid of labels on projection plots
        ax[1,1].tick_params(axis='both', which='both', 
                            labelleft=False, labelbottom=False)
        ax[0,1].tick_params(axis='both', which='both', 
                            labelleft=False, labelbottom=False)
        ax[1,0].tick_params(axis='both', which='both', 
                            labelleft=False, labelbottom=False)


        if show_ticks:
            # place ticks along shared axes as white to be visible
            
            # place ticks on left and top
            ax[1,1].tick_params(axis='both', which='both', 
                                left=True, top=True, bottom=False, right=False, 
                                color='white')
            
            # place ticks on left and bottom
            ax[0,1].tick_params(axis='both', which='both', 
                                left=True, top=False, bottom=True, right=False,  
                                color='white')
            
            # place ticks on right and top
            ax[1,0].tick_params(axis='both', which='both', 
                                left=False, top=True, bottom=False, right=True, 
                                color='white')
        else:
            # get rid of ticks
            ax[1,1].tick_params(axis='both', which='both', left=False, 
                                bottom=False)
            ax[0,1].tick_params(axis='both', which='both', left=False, 
                                bottom=False)
            ax[1,0].tick_params(axis='both', which='both', left=False, 
                                bottom=False)

        # set annotating text to label each projection
        rhox_str = r'$\rho_x$'
        rhoy_str = r'$\rho_y$'
        rhoz_str = r'$\rho_z$'
        rhox_str_coords = self.ymin + 0.05*self.ydim, self.zmin + 0.9*self.zdim
        rhoy_str_coords = self.zmin + 0.95*self.zdim, self.xmin + 0.1*self.xdim
        rhoz_str_coords = self.ymin + 0.05*self.ydim, self.xmin + 0.1*self.xdim
        ax[0,1].annotate(rhox_str, xy=rhox_str_coords, color='white', fontsize=35)
        ax[1,0].annotate(rhoy_str, xy=rhoy_str_coords, color='white', fontsize=35)
        ax[1,1].annotate(rhoz_str, xy=rhoz_str_coords, color='white', fontsize=35)
        
        # add time info if wanted
        if time is not None:
            nsnap, age, redshift = time
            time_title_str = f'snap: {nsnap:.0f} - '
            time_title_str += f'age: {age:.2f} kyr - '
            time_title_str += f'z = {redshift:.2f}'
            ax[0,1].set_title(time_title_str, fontsize=18)

        # add phase space colorbar
        fig.subplots_adjust(left=0.1, top=0.9)
        # cbar_ax = fig.add_axes([0.15, 0.85, 0.28, 0.015])
        cbar_ax = fig.add_axes([0.15, 0.52, 0.28, 0.015])
        fig.colorbar(im0, cax=cbar_ax, orientation="horizontal")
        cbar_ax.xaxis.set_ticks_position('top')

        # add shared colorbar to right
        fig.subplots_adjust(right=0.9, bottom=0.05)
        cbar_ax = fig.add_axes([1., 0.035, 0.05, 0.89])
        fig.colorbar(im3, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')

        fig.tight_layout()
        
        if fname is None:
            plt.show()
        else:
            plt.savefig(fname, dpi=self.dpi, bbox_inches = "tight")
        
        plt.close()
