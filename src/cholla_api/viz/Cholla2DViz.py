'''
Plotting for a 2D Cholla Simulation run

The general idea is that each simulation run AND each plotting value is going
    to have different formats.
    
    
The important variable that can be changed is plt_fmt. This is meant to be set
    within each plotting value function with its own plotting formatting
    function (plot_pressure has its own pressure_format function).
    
    The idea is that if you see something you'd like to change after making a plot,
    then you create a new key for plt_fmt with different values for whatever run
    or values you're plotting. Within the plotting function, you change those
    values. Within __init__, you place a default value.

    For example, if I wanted to change the title fontsize in density, I would 
    add "title_fs" as a key to plt_fmt, and assign specific font sizes to different
    runs within the density_format function.

Could be used standalone, but this is designed to be encapsulated within ChollaViz

In the future, I would like to make a completely modular plotting 2D function such
    that you call Cholla2DViz.value() which takes in the plotting formatting for
    that value. Then it calls a modular 2D plotting function so that all the plots
    look similar to one another, and is consistent. Current placeholder function
    name is plot_val_2D
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


class Cholla2DViz:
    '''
    Assumption that the 2 dimensions used in the simulation are x and y
        results in using the plotting data as data[:,:,0]
    '''
    
    def __init__(self, chollasnap, test_name=""):
        self.ch_snap = chollasnap
        self.test_name = test_name
        
        self.def_clb_fmt = reg3_fmt
        self.def_figsize = (10,8)
        self.def_compare_figsize = (12,6)
        

    def pressure_fmt(self, compare=False):
        '''
        set the formatting for the pressure plotting function / plot_pressure_2D
        '''
        plt_fmt = {}
        
        if compare:
            plt_fmt["fig_size"]=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        else:
            plt_fmt["fig_size"]=self.def_figsize
        
        plt_fmt["value_key"] = "pressure"
        plt_fmt["title"] = "Pressure"
        
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        
        return plt_fmt
    
    
    def density_fmt(self, compare=False):
        plt_fmt = {}
        
        if compare:
            plt_fmt["fig_size"]=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        else:
            plt_fmt["fig_size"]=self.def_figsize
        
        plt_fmt["value_key"] = "density"
        plt_fmt["title"] = "Density"
        
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    
    def velx_fmt(self, compare=False):
        plt_fmt = {}
        
        if compare:
            plt_fmt["fig_size"]=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        else:
            plt_fmt["fig_size"]=self.def_figsize
        
        plt_fmt["value_key"] = "vel_x"
        plt_fmt["title"] = "Velocity (x)"
        
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        
        return plt_fmt
    
    def vely_fmt(self, compare=False):
        plt_fmt = {}
        
        if compare:
            plt_fmt["fig_size"]=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        else:
            plt_fmt["fig_size"]=self.def_figsize
        
        plt_fmt["value_key"] = "vel_y"
        plt_fmt["title"] = "Velocity (y)"
        
        if (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        
        return plt_fmt
    
    def vel_fmt(self, compare=False):
        plt_fmt = {}
        
        if compare:
            plt_fmt["fig_size"]=self.def_compare_figsize
            if (self.test_name == "Rayleigh_Taylor"):
                plt_fmt["fig_size"]=(8,6)
        else:
            plt_fmt["fig_size"]=self.def_figsize
        
        plt_fmt["value_key"] = "vel_mag"
        plt_fmt["title"] = "Velocity"
        
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def pressure(self, plt_kwargs):
        plt_fmt = self.pressure_fmt()
        plot_value_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def density(self, plt_kwargs):
        plt_fmt = self.density_fmt()
        plot_value_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def velocity_x(self, plt_kwargs):
        plt_fmt = self.velx_fmt()
        plot_value_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def velocity_y(self, plt_kwargs):
        plt_fmt = self.vely_fmt()
        plot_value_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def velocity(self, plt_kwargs):
        plt_fmt = self.vel_fmt()
        plot_value_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def density_compare(self, ch_snap2, plt_kwargs):
        '''
        comparing density at different snapshot times
        '''
        plt_fmt = self.density_fmt(compare=True)
        plot_value_compare_2D(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def pressure_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.pressure_fmt(compare=True)
        plot_value_compare_2D(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocityx_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.velx_fmt(compare=True)
        plot_value_compare_2D(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocityy_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.vely_fmt(compare=True)
        plot_value_compare_2D(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocity_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.vel_fmt(compare=True)
        plot_value_compare_2D(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    