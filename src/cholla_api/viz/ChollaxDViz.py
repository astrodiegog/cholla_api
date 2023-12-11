'''
Helper plotting function that handles some of the similar things that
    come up with visualizations of all dimensions
'''

from cholla_api.viz.Cholla1DViz import Cholla1DVizFmt
from cholla_api.viz.Cholla2DViz import Cholla2DVizFmt


class ChollaxDViz:
    
    def __init__(self, chollasnap, num_dims, test_name=""):
        self.ch_snap = chollasnap
        self.num_dims = num_dims
        self.test_name = test_name
        
        if (self.num_dims == 1):
            ch_1DVizFmt = Cholla1DVizFmt(test_name=self.test_name)
            self.plot_value = ch_1DVizFmt.plot_value
            self.plot_value_compare = ch_1DVizFmt.plot_value_compare
            
            self.pressure_fmt = ch_1DVizFmt.pressure_fmt
            self.density_fmt = ch_1DVizFmt.density_fmt
            self.vel_fmt = ch_1DVizFmt.vel_fmt
            
            self.pressure_compare_fmt = ch_1DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_1DVizFmt.density_fmt
            self.vel_compare_fmt = ch_1DVizFmt.vel_fmt
            
        elif (self.num_dims == 2):
            ch_2DVizFmt = Cholla2DVizFmt(test_name=self.test_name)
            
            self.plot_value = ch_2DVizFmt.plot_value
            self.plot_value_compare = ch_2DVizFmt.plot_value_compare
            
            self.pressure_fmt = ch_2DVizFmt.pressure_fmt
            self.density_fmt = ch_2DVizFmt.density_fmt
            self.velx_fmt = ch_2DVizFmt.velx_fmt
            self.vely_fmt = ch_2DVizFmt.vely_fmt
            self.vel_fmt = ch_2DVizFmt.vel_fmt
            
            self.pressure_compare_fmt = ch_2DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_2DVizFmt.density_fmt
            self.velx_compare_fmt = ch_2DVizFmt.velx_fmt
            self.vely_compare_fmt = ch_2DVizFmt.vely_fmt
            self.vel_compare_fmt = ch_2DVizFmt.vel_fmt
        
            
    def pressure(self, plt_kwargs):
        plt_fmt = self.pressure_fmt()
        self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
    
    def density(self, plt_kwargs):
        plt_fmt = self.density_fmt()
        self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def velocity(self, plt_kwargs):
        plt_fmt = self.vel_fmt()
        self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        
    def velocity_x(self, plt_kwargs):
        plt_fmt = self.velx_fmt()
        self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
    
    def velocity_y(self, plt_kwargs):
        plt_fmt = self.vely_fmt()
        self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
    
    
    def density_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.density_compare_fmt(compare=True)
        self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def pressure_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.pressure_compare_fmt(compare=True)
        self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocity_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.vel_compare_fmt(compare=True)
        self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocityx_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.velx_compare_fmt(compare=True)
        self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
    
    def velocityy_compare(self, ch_snap2, plt_kwargs):
        plt_fmt = self.vely_compare_fmt(compare=True)
        self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, 
                              ch_snap2.data, ch_snap2.head, plt_fmt, plt_kwargs)
        



