'''
Helper plotting function that handles some of the similar things that
    come up with visualizations of all dimensions
'''

from cholla_api.viz.Cholla1DViz import Cholla1DVizFmt
from cholla_api.viz.Cholla2DViz import Cholla2DVizFmt
from cholla_api.viz.Cholla3DViz import Cholla3DVizFmt

class ChollaxDViz:
    '''
    Helper class that sets plotting functions wrt their dimensions
    '''
    def __init__(self, chollasnap, num_dims, test_name="", reduce_dim_fn=None):
        self.ch_snap = chollasnap
        self.num_dims = num_dims
        self.test_name = test_name
        self.reduce_dim_fn = reduce_dim_fn
        
        # create the correct dimension visualization format, assign plotting functions
        if (self.num_dims == 1):
            ch_1DVizFmt = Cholla1DVizFmt(test_name=self.test_name)
            self.plot_value = ch_1DVizFmt.plot_value
            self.plot_value_compare = ch_1DVizFmt.plot_value_compare
            self.plot_value_reduce = None # nonsense to reduce 1D data
            
            self.pressure_fmt = ch_1DVizFmt.pressure_fmt
            self.density_fmt = ch_1DVizFmt.density_fmt
            self.vel_fmt = ch_1DVizFmt.vel_fmt
            
            self.pressure_compare_fmt = ch_1DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_1DVizFmt.density_fmt
            self.vel_compare_fmt = ch_1DVizFmt.vel_fmt
            
        elif (self.num_dims == 2):
            ch_2DVizFmt = Cholla2DVizFmt(test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
            
            self.plot_value = ch_2DVizFmt.plot_value
            self.plot_value_compare = ch_2DVizFmt.plot_value_compare
            self.plot_value_reduce = ch_2DVizFmt.plot_value_reduce
            
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
        
        elif (self.num_dims == 3):
            ch_3DVizFmt = Cholla3DVizFmt(test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)

            self.plot_value = ch_3DVizFmt.plot_value
            self.plot_value_compare = ch_3DVizFmt.plot_value_compare
            self.plot_value_reduce = ch_3DVizFmt.plot_value_reduce

            self.pressure_fmt = ch_3DVizFmt.pressure_fmt
            self.density_fmt = ch_3DVizFmt.density_fmt
            self.velx_fmt = ch_3DVizFmt.velx_fmt
            self.vely_fmt = ch_3DVizFmt.vely_fmt
            self.vel_fmt = ch_3DVizFmt.vel_fmt

            self.pressure_compare_fmt = ch_3DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_3DVizFmt.density_fmt
            self.velx_compare_fmt = ch_3DVizFmt.velx_fmt
            self.vely_compare_fmt = ch_3DVizFmt.vely_fmt
            self.vel_compare_fmt = ch_3DVizFmt.vel_fmt
            
    def make_plot(self, plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None):
        '''
        given the plotting format required for the data, make the plotting call
        '''
        if reduce_dim:
            self.plot_value_reduce(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)
        elif compare_snap2:
            self.plot_value_compare(self.ch_snap.data, self.ch_snap.head, compare_snap2.data,
                                    compare_snap2.head, plt_fmt, plt_kwargs)
        else:
            self.plot_value(self.ch_snap.data, self.ch_snap.head, plt_fmt, plt_kwargs)

    def pressure(self, plt_kwargs, reduce_dim=None):
        '''
        sets pressure plotting formats and creates plot
        '''
        plt_fmt = self.pressure_fmt()
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None)
    
    def density(self, plt_kwargs, reduce_dim=None):
        '''
        sets density plotting formats and creates plot
        '''
        plt_fmt = self.density_fmt()
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None)
        
    def velocity(self, plt_kwargs, reduce_dim=None):
        '''
        sets velocity plotting formats and creates plot
        '''
        plt_fmt = self.vel_fmt()
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None)
 
    def velocity_x(self, plt_kwargs, reduce_dim=None):
        '''
        sets x velocity plotting formats and creates plot
        '''
        plt_fmt = self.velx_fmt()
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None)
    
    def velocity_y(self, plt_kwargs, reduce_dim=None):
        '''
        sets y velocity plotting formats and creates plot
        '''
        plt_fmt = self.vely_fmt()
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=None)
    
    
    def density_compare(self, ch_snap2, plt_kwargs, reduce_dim=None):
        '''
        sets comparing density plotting formats and creates plot
        '''
        plt_fmt = self.density_compare_fmt(compare=True)
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=ch_snap2)
    
    def pressure_compare(self, ch_snap2, plt_kwargs, reduce_dim=None):
        '''
        sets comparing pressure plotting formats and creates plot
        '''
        plt_fmt = self.pressure_compare_fmt(compare=True)
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=ch_snap2)
 
    def velocity_compare(self, ch_snap2, plt_kwargs, reduce_dim=None):
        '''
        sets comparing velocity plotting formats and creates plot
        '''
        plt_fmt = self.vel_compare_fmt(compare=True)
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=ch_snap2)
 
    def velocityx_compare(self, ch_snap2, plt_kwargs, reduce_dim=None):
        '''
        sets comparing x velocity plotting formats and creates plot
        '''
        plt_fmt = self.velx_compare_fmt(compare=True)
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=ch_snap2)
 
    def velocityy_compare(self, ch_snap2, plt_kwargs, reduce_dim=None):
        '''
        sets comparing y velocity plotting formats and creates plot
        '''
        plt_fmt = self.vely_compare_fmt(compare=True)
        self.make_plot(plt_fmt, plt_kwargs, reduce_dim, compare_snap2=ch_snap2) 



