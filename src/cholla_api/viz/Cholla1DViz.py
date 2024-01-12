'''
Plotting for a 1D Cholla Simulation run
'''

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

from cholla_api.viz.viz_format import *



def plot_value_1D(data, head, plt_fmt, plt_kwargs):
    '''
    generic plotting 1D plotting function given
    params:
        data (np array): data from cholla snapshot
        head (dictionary): header from cholla snapshot
        plt_fmt (dictionary): handles specific value formats for different tests
        plt_kwargs (dictionary): handles plot outputs (showing, saving)
    '''
    val_key = plt_fmt["value_key"]
    title_str = plt_fmt["title"]
    yval_fmt = plt_fmt["yval_fmt"]
    figsize = plt_fmt["fig_size"]
    
    value = data[val_key][:,0,0].T
    time = head["t"][0]
    time_str = f"t = {time:.3f}"
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    im = ax.plot(value)

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(yval_fmt))

    ax.set_xlim(0, head["dims"][0])
    if plt_fmt.get("ylims"):
        ax.set_ylim(plt_fmt["ylims"])
    
    ax.set_xlabel("X (cells)")
    ax.set_ylabel(title_str)
    ax.set_title(title_str)
    ax.text(0.8, 1.02, time_str, fontsize=26, transform=ax.transAxes)
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()


def plot_value_compare_1D(data1, head1, data2, head2, plt_fmt, plt_kwargs):
    '''
    generic plotting 1D plotting comparison function
    params:
        data1/data2 (np array): data from cholla snapshot
        head1/head2 (dictionary): header from cholla snapshot
        plt_fmt (dictionary): handles specific value formats for different tests
        plt_kwargs (dictionary): handles plot outputs (showing, saving)
    '''
    val_key = plt_fmt["value_key"]
    y_label = plt_fmt["y_label"]
    yval_fmt = plt_fmt["yval_fmt"]
    figsize = plt_fmt["fig_size"]
    
    value1 = data1[val_key][:,0,0].T
    value2 = data2[val_key][:,0,0].T
    values = [value1, value2]
    
    t1 = head1["t"][0]
    t2 = head2["t"][0]
    all_t = [t1,t2]
    
    xlims = (0, head1["dims"][0])
    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

    for i in range(2):
        ax[i].plot(values[i])

        ax[i].set_xlim(xlims)

        ax[i].set_xlabel("X (cells)")
        ax[i].set_title(f"t = {all_t[i]:.3f}", fontsize=26)

        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(yval_fmt))

    ax[0].set_ylabel(y_label)
    
    plt.tight_layout()
    if plt_kwargs.get("imgfout"):
        plt.savefig(plt_kwargs["imgfout"])
    if plt_kwargs.get("show"):
        plt.show()
    plt.close()



class Cholla1DVizFmt:
    '''
    Assumption that the 1 dimension used in the simulation is x 
        results in using the plotting data as data[:,0,0]
    '''
    

    def __init__(self, test_name=""):

        self.test_name = test_name
        
        self.def_yval_fmt = reg1_fmt
        self.def_figsize = (10,8)
        self.def_compare_figsize = (12,6)
        
        self.plot_value = plot_value_1D
        self.plot_value_compare = plot_value_compare_1D
        

    def pressure_fmt(self, compare=False):
        '''
        set the formatting for the pressure plotting function
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "pressure"
        
        # figsize & (y label OR title)
        figsize=self.def_figsize
        if compare:
            plt_fmt["y_label"] = "Pressure"
            figsize = self.def_compare_figsize
        else:
            plt_fmt["title"] = "Pressure"
        plt_fmt["fig_size"] = figsize
        
        # y value
        yval_fmt = self.def_yval_fmt
        if (self.test_name == "123"):
            yval_fmt = reg5_fmt
        plt_fmt["yval_fmt"] = yval_fmt

        # y lims
        ylim_fmt = (0,1)
        if (self.test_name == "123"):
            ylim_fmt = (0,0.4)
        elif (self.test_name == "blast"):
            ylim_fmt = (0,1000)
        elif (self.test_name == "constant" or self.test_name == "square_wave"):
            ylim_fmt = (-0.1,0.1)
        elif (self.test_name == "strong_shock"):
            ylim_fmt = (0,100)
        elif (self.test_name == "two_shocks"):
            ylim_fmt = (0,2000)
        plt_fmt["ylims"] = ylim_fmt
        
        return plt_fmt
    
    def density_fmt(self, compare=False):
        '''
        set the formatting for the density plotting function
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "density"
        
        # figsize & (y label OR title)
        figsize = self.def_figsize
        if compare:
            plt_fmt["y_label"] = "Density"
            figsize = self.def_compare_figsize
        else:
            plt_fmt["title"] = "Density"
        plt_fmt["fig_size"] = figsize
        
        # y value
        yval_fmt = self.def_yval_fmt
        if (self.test_name == "123"):
            yval_fmt = reg3_fmt
        elif (self.test_name == "constant"):
            yval_fmt = reg0_fmt
        plt_fmt["yval_fmt"] = yval_fmt

        # y lims
        ylim_fmt = (0,1)
        if (self.test_name == "blast"):
            ylim_fmt = (0,8.)
        elif (self.test_name == "constant"):
            ylim_fmt = (9999, 10001)
        elif (self.test_name == "square_wave"):
            ylim_fmt = (0.8, 1.6)
        elif (self.test_name == "strong_shock"):
            ylim_fmt = (0,10)
        elif (self.test_name == "two_shocks"):
            ylim_fmt = (0,30)
        plt_fmt["ylims"] = ylim_fmt

        return plt_fmt
    
    def vel_fmt(self, compare=False):
        '''
        set the formatting for the velocity plotting function
        '''
        plt_fmt = {}
        plt_fmt["value_key"] = "vel_mag"
        
        # figsize & (y label OR title)
        figsize=self.def_figsize
        if compare:
            plt_fmt["y_label"] = "Velocity"
            figsize = self.def_compare_figsize
        else:
            plt_fmt["title"] = "Velocity"
        plt_fmt["fig_size"] = figsize
        
        # y value
        yval_fmt = self.def_yval_fmt
        if (self.test_name == "123"):
            yval_fmt = reg3_fmt
        plt_fmt["yval_fmt"] = yval_fmt
        
        # y lims
        ylim_fmt = (0,22)
        if (self.test_name == "123"):
            ylim_fmt = (0,2)
        elif (self.test_name == "constant"):
            ylim_fmt = (-1,1)
        elif (self.test_name == "sod"):
            ylim_fmt = (0,1)
        elif (self.test_name == "square_wave"):
            ylim_fmt = (0,2)
        elif (self.test_name == "strong_shock"):
            ylim_fmt = (0,5)
        elif (self.test_name == "trac_pen"):
            ylim_fmt = (118,122)
        plt_fmt["ylims"] = ylim_fmt

        return plt_fmt

