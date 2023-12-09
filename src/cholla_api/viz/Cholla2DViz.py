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


def plot_pressure_2D(data, head, plt_fmt, plt_kwargs):
    
    clb_fmt = plt_fmt["colorbar_fmt"]

    pressure = data["pressure"][:,:,0].T
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1)

    im = ax.imshow(pressure)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)


    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title("Pressure")
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()

    
def plot_density_2D(data, head, plt_fmt, plt_kwargs):
    
    clb_fmt = plt_fmt["colorbar_fmt"]

    density = data["density"][:,:,0].T
    time = head["t"][0]
    time_str = f"t = {time:.3f}"
    
    fig, ax = plt.subplots(nrows=1, ncols=1)

    im = ax.imshow(density)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title("Density")
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()
    
    
    
def plot_velx_2D(data, head, plt_fmt, plt_kwargs):

    clb_fmt = plt_fmt["colorbar_fmt"]

    vel_x = data["vel_x"][:,:,0].T
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1)

    im = ax.imshow(vel_x)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title("Velocity - x")
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()


def plot_vely_2D(data, head, plt_fmt, plt_kwargs):
    
    clb_fmt = plt_fmt["colorbar_fmt"]

    vel_y = data["vel_y"][:,:,0].T
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1)

    im = ax.imshow(vel_y)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title("Velocity - y")
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()

def plot_vel_2D(data, head, plt_fmt, plt_kwargs):
    
    clb_fmt = plt_fmt["colorbar_fmt"]
    
    vel = data["vel_mag"][:,:,0].T
    
    time = head["t"][0]
    time_str = f"t = {time:.3f}"

    fig, ax = plt.subplots(nrows=1, ncols=1)

    im = ax.imshow(vel)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im, cax=cax, orientation='vertical', format=ticker.FuncFormatter(clb_fmt))

    ax.set_xlim(0, head["dims"][0] - 1)
    ax.set_ylim(0, head["dims"][1] - 1)
    
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel("Y (cells)")
    ax.set_title("Velocity")
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()
    

    
def plot_val_2D():
    '''
    would like to make a general 2D plotting format where I can just pass in the
        data/head values and any formatting details (plt_fmt), and have a 2D plot
    '''
    return


def plot_val_1D():
    '''
    same as plot_val_2D, except looking at 1D. would also pass along some sort of
        dimension_reduction function (median along x, sum along y)
    '''
    return


class Cholla2DViz:
    '''
    Assumption that the 2 dimensions used in the simulation are x and y
    '''
    
    def __init__(self, chollasnap, test_name=""):
        self.ch_snap = chollasnap
        self.test_name = test_name
        
        self.def_clb_fmt = reg3_fmt
        

    def pressure_fmt(self):
        '''
        set the formatting for the pressure plotting function / plot_pressure_2D
        '''
        plt_fmt = {}
        
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def pressure(self, plt_kwargs):
        plt_fmt = self.pressure_fmt()
        plot_pressure_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
    

    def density_fmt(self):
        '''
        set the formatting for the pressure plotting function / plot_density_2D
        '''
        plt_fmt = {}
        
        if (self.test_name == "KH_resind"):
            clb_fmt = reg1_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = reg5_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def density(self, plt_kwargs):
        plt_fmt = self.density_fmt()
        plot_density_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
        
    def velx_fmt(self):
        '''
        set the formatting for the pressure plotting function / plot_velx_2D
        '''
        plt_fmt = {}
        
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def velocity_x(self, plt_kwargs):
        plt_fmt = self.velx_fmt()
        plot_velx_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
        
    def vely_fmt(self):
        '''
        set the formatting for the pressure plotting function / plot_vely_2D
        '''
        plt_fmt = {}
        
        if (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def velocity_y(self, plt_kwargs):
        plt_fmt = self.vely_fmt()
        plot_vely_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    
    def vel_fmt(self):
        '''
        set the formatting for the pressure plotting function / plot_vely_2D
        '''
        plt_fmt = {}
        
        if (self.test_name == "sod"):
            clb_fmt = scinot2_fmt
        elif (self.test_name == "sound_wave"):
            clb_fmt = scinot2_fmt
        else:
            clb_fmt = self.def_clb_fmt
        plt_fmt["colorbar_fmt"] = clb_fmt
        
        return plt_fmt
    
    def velocity(self, plt_kwargs):
        plt_fmt = self.vely_fmt()
        plot_vel_2D(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    