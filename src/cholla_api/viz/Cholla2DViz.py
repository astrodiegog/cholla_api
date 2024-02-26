'''
Plotting for a 2D Cholla Simulation run
'''

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

from cholla_api.viz.viz_format import *



def plot_value_2D(data, head, plt_fmt, plt_kwargs):
    '''
    generic plotting 2D plotting function given
    params:
        data (np array): data from cholla snapshot
        head (dictionary): header from cholla snapshot
        plt_fmt (dictionary): handles specific value formats for different tests
        plt_kwargs (dictionary): handles plot outputs (showing, saving)
    '''    
    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    clb_fmt = plt_fmt["value_fmt"]
    figsize = plt_fmt["fig_size"]
    
    value = data[val_key][:,:,0].T
    
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    if (plt_fmt.get("value_lims")):
        vmin, vmax = plt_fmt["value_lims"]
    else:
        vmin, vmax = np.min(value), np.max(value)

    im = ax.imshow(value, vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title(title_str)
    plt.tight_layout()
    ax.text(0.8, 1.02, time_str, fontsize=26, transform=ax.transAxes)
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()

def plot_value_compare_2D(data1, head1, data2, head2, plt_fmt, plt_kwargs):
    '''
    generic 2D plotting comparison function given
    params:
        data1/data2 (np array): data from cholla snapshot
        head1/head2 (dictionary): header from cholla snapshot
        plt_fmt (dictionary): handles specific value formats for different tests
        plt_kwargs (dictionary): handles plot outputs (showing, saving)
    '''
    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    clb_fmt = plt_fmt["value_fmt"]
    figsize = plt_fmt["fig_size"]
    
    value1 = data1[val_key][:,:,0].T
    value2 = data2[val_key][:,:,0].T
    values = [value1, value2]
    
    t1 = head1["t"][0]
    t2 = head2["t"][0]
    all_t = [t1,t2]
    
    xlims = (0, head1["dims"][0] - 1)
    ylims = (0, head1["dims"][1] - 1)
    
    gridspec = {'width_ratios': [1, 1, 0.08]}
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=figsize, gridspec_kw=gridspec)
    if (plt_fmt.get("value_lims")):
        vmin, vmax = plt_fmt["value_lims"]
    else:
        vmin, vmax = np.min(values), np.max(values)


    for i in range(2):
        im = ax[i].imshow(values[i], vmin=vmin, vmax=vmax)

        ax[i].set_xlim(xlims)
        ax[i].set_ylim(ylims)

        ax[i].set_xlabel("X (cells)")
        ax[i].set_title(f"t = {all_t[i]:.4f}", fontsize=28)

        
    cax = ax[2]
    fig.colorbar(im, cax=cax, orientation='vertical', 
                 format=ticker.FuncFormatter(clb_fmt))

    
    ax[1].get_yaxis().set_ticks([])
    ax[0].set_ylabel("Y (cells)")
    
    fig.suptitle(title_str)
    fig.tight_layout()
    fig.subplots_adjust(top=0.8, wspace=0.05)
    
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()
    
    
def def_reduce_fn(data):
    '''
    default dimension reduction function. median along axis=1 for a 2D data
    '''
    return np.median(data, axis=1)


def plot_value_2D_1D(data, head, plt_fmt, plt_kwargs, red_fn=def_reduce_fn):
    '''
    same as plot_value_2D, but accepts a dimension reduction function
        for example, median along x, sum along y
        default would probably be just taking the median
    '''

    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    yval_fmt = plt_fmt["value_fmt"]    
    figsize = plt_fmt["fig_size"]

    value = data[val_key][:,:,0].T
    value_reduced = red_fn(value) # make sure default_reduce works

    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    
    im = ax.plot(value_reduced)

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(yval_fmt))

    ax.set_xlim(0, head["dims"][0])
    if (plt_fmt.get("value_lims")):
        vmin, vmax = plt_fmt["value_lims"]
    else:
        vmin, vmax = np.min(value), np.max(value)
    ax.set_ylim(vmin, vmax)

    ax.set_xlabel("X (cells)")
    if plt_fmt.get("x_label"):
        ax.set_xlabel(plt_fmt["x_label"])

    ax.set_ylabel(title_str)
    ax.set_title(title_str)
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()




class Cholla2DVizFmt:
    '''
    Assumption that the 2 dimensions used in the simulation are x and y
        results in using the plotting data as data[:,:,0]
    '''
    
    def __init__(self, test_name="", reduce_dim_fn=None):
        self.test_name = test_name
        self.reduce_dim_fn = reduce_dim_fn
        
        self.def_val_fmt = reg3_fmt
        self.def_figsize = (10,8)
        self.def_compare_figsize = (12,6)

        self.def_reduce_fn = def_reduce_fn
        # make sure reduce fn is a function if it is to be used
        if hasattr(self.reduce_dim_fn, '__call__'):
            self.def_reduce_fn = self.reduce_dim_fn
        
        self.plot_value = plot_value_2D
        self.plot_value_compare = plot_value_compare_2D
        self.plot_value_reduce = plot_value_2D_1D
        

    def pressure_fmt(self, compare=False):
        '''
        set the plotting function formatting for pressure
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "pressure"
        plt_fmt["title"] = "Pressure"
        
        # figsize
        fig_size=self.def_figsize
        if compare:
            fig_size=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                fig_size=(8,6)
        plt_fmt["fig_size"] = fig_size
        
        # colorbar/y-value for dimension reduction
        val_fmt = self.def_val_fmt
        if (self.test_name == "KH_resind"):
            val_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            val_fmt = reg5_fmt
        plt_fmt["value_fmt"] = val_fmt

        # value limits
        val_lims = (1.5, 3.)
        if (self.test_name == "implosion"):
            val_lims = (0.5, 1.5)
        elif (self.test_name == "sod"):
            val_lims = (0., 1.)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (0.5, 1.0)
        elif (self.test_name == "sound_wave"):
            val_lims = (0.5999,0.6001)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt
    
    
    def density_fmt(self, compare=False):
        '''
        set the plotting function formattig for density
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "density"
        plt_fmt["title"] = "Density"
        
        # figsize
        fig_size=self.def_figsize
        if compare:
            fig_size=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                fig_size=(8,6)
        plt_fmt["fig_size"]=fig_size
        
        # colorbar/y-value
        value_fmt = self.def_val_fmt
        if (self.test_name == "KH_resind"):
            value_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            value_fmt = reg5_fmt
        plt_fmt["value_fmt"] = value_fmt
        
        # value limits
        val_lims = (0.,1.)
        if (self.test_name == "KH_disc"):
            val_lims = (1., 2.)
        elif (self.test_name == "KH_resind"):
            val_lims = (0., 100.)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (1., 2.)
        elif (self.test_name == "sound_wave"):
            val_lims = (0.9999, 1.0001)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt
    
    
    def velx_fmt(self, compare=False):
        '''
        set the plotting function formattig for x velocity
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "vel_x"
        plt_fmt["title"] = "Velocity (x)"
        
        # figsize
        figsize=self.def_figsize
        if compare:
            figsize=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                figsize=(8,6)   
        plt_fmt["fig_size"] = figsize
        
        # colorbar/y-value
        value_fmt = self.def_val_fmt
        if (self.test_name == "sound_wave"):
            value_fmt = scinot2_fmt
        plt_fmt["value_fmt"] = value_fmt
        
        # value limits
        val_lims = (0.,1.)
        if (self.test_name == "implosion"):
            val_lims = (-0.15, 0.15)
        elif (self.test_name == "KH_disc"):
            val_lims = (-1.,1.)
        elif (self.test_name == "KH_resind"):
            val_lims = (-0.5, 0.5)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (-0.002, 0.002)
        elif (self.test_name == "sound_wave"):
            val_lims = (-6e-5,6e-5)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt
    
    def vely_fmt(self, compare=False):
        '''
        set the plotting function formatting for y velocity
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "vel_y"
        plt_fmt["title"] = "Velocity (y)"
        
        # figsize
        fig_size = self.def_figsize
        if compare:
            fig_size=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        plt_fmt["fig_size"]=fig_size
        
        # colorbar/y-value
        value_fmt = self.def_val_fmt
        if (self.test_name == "sound_wave"):
            value_fmt = scinot2_fmt
        plt_fmt["value_fmt"] = value_fmt
        
        # value limits
        val_lims = (-1.,1.)
        if (self.test_name == "implosion"):
            val_lims = (-0.15, 0.15)
        elif (self.test_name == "KH_resind"):
            val_lims = (-0.5, 0.5)
        elif (self.test_name == "sod"):
            val_lims = (-0.1, 0.1)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (-0.015, 0.015)
        elif (self.test_name == "sound_wave"):
            val_lims = (-6e-5,6e-5)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt
    
    def vel_fmt(self, compare=False):
        '''
        set the plotting function formatting for velocity
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "vel_mag"
        plt_fmt["title"] = "Velocity"
        
        # figsize
        fig_size=self.def_figsize
        if compare:
            fig_size=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        plt_fmt["fig_size"]=fig_size
        
        # colorbar
        value_fmt = self.def_val_fmt
        if (self.test_name == "sound_wave"):
            value_fmt = scinot2_fmt
        plt_fmt["value_fmt"] = value_fmt
        
        # value limits
        val_lims = (0.,1.)
        if (self.test_name == "implosion"):
            val_lims = (-0.15, 0.15)
        elif (self.test_name == "KH_resind"):
            val_lims = (0., 0.5)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (0., 0.02)
        elif (self.test_name == "sound_wave"):
            val_lims = (0, 2e-4)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt


    def inten_fmt(self, compare=False):
        '''
        set the formatting for the internal energy
        '''
        
        plt_fmt = {}
        plt_fmt["value_key"] = "int_energy"
        plt_fmt["title"] = "Int. Energy"
        
        # figsize
        fig_size=self.def_figsize
        if compare:
            fig_size=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        plt_fmt["fig_size"]=fig_size
        
        # colorbar
        value_fmt = self.def_val_fmt
        if (self.test_name == "sound_wave"):
            value_fmt = scinot3_fmt
            
        plt_fmt["value_fmt"] = value_fmt

        # value limits
        val_lims = (0,4.)
        if (self.test_name == "implosion"):
            val_lims = (2., 5.)
        elif (self.test_name == "KH_disc"):
            val_lims = (2., 4.)
        elif (self.test_name == "sod"):
            val_lims = (1.5, 3.5)
        elif (self.test_name == "Rayleigh_Taylor"):
            val_lims = (0.5, 2.)
        elif (self.test_name == "sound_wave"):
            val_lims = (8e-1, 1.)
        plt_fmt["value_lims"] = val_lims

        return plt_fmt
