'''
Plotting for a 2D Cholla Simulation run
'''

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

from cholla_api.viz.viz_format import *



def plot_value_2D(data, head, plt_fmt, plt_kwargs):
    
    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    clb_fmt = plt_fmt["colorbar_fmt"]
    figsize = plt_fmt["fig_size"]
    
    value = data[val_key][:,:,0].T
    
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    im = ax.imshow(value)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title(title_str)
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()

def plot_value_compare_2D(data1, head1, data2, head2, plt_fmt, plt_kwargs):
    
    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    clb_fmt = plt_fmt["colorbar_fmt"]
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
    
    for i in range(2):
        im = ax[i].imshow(values[i])

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
    
    

def plot_value_2D_1D():
    '''
    same as plot_value_2D, but accepts a dimension reduction function
        for example, median along x, sum along y
        default would probably be just taking the median
    '''
    return


class Cholla2DVizFmt:
    '''
    Assumption that the 2 dimensions used in the simulation are x and y
        results in using the plotting data as data[:,:,0]
    '''
    
    def __init__(self, test_name=""):
        self.test_name = test_name
        
        self.def_clb_fmt = reg3_fmt
        self.def_figsize = (10,8)
        self.def_compare_figsize = (12,6)
        
        self.plot_value = plot_value_2D
        self.plot_value_compare = plot_value_compare_2D
        

    def pressure_fmt(self, compare=False):
        '''
        set the formatting for the pressure plotting function / plot_pressure_2D
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
        
        # colorbar
        clb_fmt = self.def_clb_fmt
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    
    def density_fmt(self, compare=False):
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
        
        # colorbar
        clb_fmt = self.def_clb_fmt
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    
    def velx_fmt(self, compare=False):
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
        
        # colorbar
        clb_fmt = self.def_clb_fmt
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def vely_fmt(self, compare=False):
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
        
        # colorbar
        clb_fmt = self.def_clb_fmt
        if (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def vel_fmt(self, compare=False):
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
        clb_fmt = self.def_clb_fmt
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt

    