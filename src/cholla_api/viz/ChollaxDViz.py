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
            self.inten_fmt = ch_1DVizFmt.inten_fmt
            
            self.pressure_compare_fmt = ch_1DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_1DVizFmt.density_fmt
            self.vel_compare_fmt = ch_1DVizFmt.vel_fmt
            self.inten_compare_fmt = ch_1DVizFmt.inten_fmt
            
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
            self.inten_fmt = ch_2DVizFmt.inten_fmt
            
            self.pressure_compare_fmt = ch_2DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_2DVizFmt.density_fmt
            self.velx_compare_fmt = ch_2DVizFmt.velx_fmt
            self.vely_compare_fmt = ch_2DVizFmt.vely_fmt
            self.vel_compare_fmt = ch_2DVizFmt.vel_fmt
            self.inten_compare_fmt = ch_2DVizFmt.inten_fmt
        
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
            self.inten_fmt = ch_3DVizFmt.inten_fmt

            self.pressure_compare_fmt = ch_3DVizFmt.pressure_fmt
            self.density_compare_fmt = ch_3DVizFmt.density_fmt
            self.velx_compare_fmt = ch_3DVizFmt.velx_fmt
            self.vely_compare_fmt = ch_3DVizFmt.vely_fmt
            self.vel_compare_fmt = ch_3DVizFmt.vel_fmt
            self.inten_compare_fmt = ch_3DVizFmt.inten_fmt
    
    def make_plot(self, dataftype, plt_fmt, plt_kwargs, compare_snap2=None):
        '''
        given the plotting format required for the data, make the plotting call
        '''
        cholla_data = None
        cholla_data2 = None
        if dataftype == "hydro":
            cholla_data = self.ch_snap.data_hydro
            if compare_snap2:
                cholla_data2 = compare_snap2.data_hydro
        elif dataftype == "parts":
            cholla_data = self.ch_snap.data_parts
            if compare_snap2:
                cholla_data2 = compare_snap2.data_parts
        elif dataftype == "gravy":
            cholla_data = self.ch_snap.data_gravy
            if compare_snap2:
                cholla_data2 = compare_snap2.data_gravy
        else:
            # handle error where file is alive or ftype not specified
            return -1
        
        if self.reduce_dim_fn:
            self.plot_value_reduce(cholla_data, self.ch_snap.head, plt_fmt, plt_kwargs, red_fn=self.reduce_dim_fn)
        elif compare_snap2:
            self.plot_value_compare(cholla_data, self.ch_snap.head, cholla_data2,
                                    compare_snap2.head, plt_fmt, plt_kwargs)
        else:
            self.plot_value(cholla_data, self.ch_snap.head, plt_fmt, plt_kwargs)
    
    def make_hydroplot(self, plt_fmt, plt_kwargs, compare_snap2=None):
        self.make_plot(dataftype="hydro", plt_fmt=plt_fmt, plt_kwargs=plt_kwargs,
                       compare_snap2=compare_snap2)
        
    def make_partsplot(self, plt_fmt, plt_kwargs, compare_snap2=None):
        self.make_plot(dataftype="parts", plt_fmt=plt_fmt, plt_kwargs=plt_kwargs,
                       compare_snap2=compare_snap2)
        
    def make_gravyplot(self, plt_fmt, plt_kwargs, compare_snap2=None):
        self.make_plot(dataftype="gravy", plt_fmt=plt_fmt, plt_kwargs=plt_kwargs,
                       compare_snap2=compare_snap2)

    def pressure(self, plt_kwargs):
        '''
        sets pressure plotting formats and creates plot
        '''
        plt_fmt = self.pressure_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
    
    def density(self, plt_kwargs):
        '''
        sets density plotting formats and creates plot
        '''
        plt_fmt = self.density_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
        
    def int_energy(self, plt_kwargs):
        '''
        sets internal energy plotting formats and creates plot
        '''
        plt_fmt = self.inten_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
        
    def velocity(self, plt_kwargs):
        '''
        sets velocity plotting formats and creates plot
        '''
        plt_fmt = self.vel_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
 
    def velocity_x(self, plt_kwargs):
        '''
        sets x velocity plotting formats and creates plot
        '''
        plt_fmt = self.velx_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
    
    def velocity_y(self, plt_kwargs):
        '''
        sets y velocity plotting formats and creates plot
        '''
        plt_fmt = self.vely_fmt()
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=None)
    
    
    def density_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing density plotting formats and creates plot
        '''
        plt_fmt = self.density_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2)
    
    def pressure_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing pressure plotting formats and creates plot
        '''
        plt_fmt = self.pressure_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2)
        
    def internal_energy_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing internal energy plotting formats and creates plot
        '''
        plt_fmt = self.inten_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2)
 
    def velocity_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing velocity plotting formats and creates plot
        '''
        plt_fmt = self.vel_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2)
 
    def velocityx_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing x velocity plotting formats and creates plot
        '''
        plt_fmt = self.velx_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2)
 
    def velocityy_compare(self, ch_snap2, plt_kwargs):
        '''
        sets comparing y velocity plotting formats and creates plot
        '''
        plt_fmt = self.vely_compare_fmt(compare=True)
        self.make_hydroplot(plt_fmt, plt_kwargs, compare_snap2=ch_snap2) 



