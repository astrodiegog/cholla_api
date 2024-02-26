'''
Handles all visualization API with user. Should be simple
'''


'''
when catching the keyerror for each value, should I be raising
    KeyError or just returning a None?
'''

import numpy as np

from cholla_api.viz.ChollaxDViz import ChollaxDViz


class ChollaViz:
    '''
    Class that handles any visualization given a Cholla snapshot.
    '''
    
    def __init__(self, chollasnap, test_name="", plt_kwargs={"show": True, "save": False}, reduce_dim_fn=None):
        self.ch_snap = chollasnap
        self.test_name = test_name
        self.plt_kwargs = plt_kwargs
        self.reduce_dim_fn = reduce_dim_fn
        self.num_D = int(np.sum(self.ch_snap.head['dims'] != 1.))
        
        
        # create the helper xDViz class that is flexible to number of dimensions
        if self.num_D == 1:
            self.chollaxDViz = ChollaxDViz(self.ch_snap, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
            
        elif self.num_D == 2:
            self.chollaxDViz = ChollaxDViz(self.ch_snap, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
            
        elif self.num_D == 3:
            self.chollaxDViz = ChollaxDViz(self.ch_snap, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
        
        self.set_value_keys()
        self.set_plot_keys()
        
    
    def set_plot_keys(self):
        '''
        helper function to find appropriate plotting method to call
        '''
        self.plot_keys = {}
        self.plot_keys[self.pressure_key] = self.pressure
        self.plot_keys[self.density_key] = self.density
        self.plot_keys[self.vel_str] = self.velocity
        self.plot_keys[self.inten_key] = self.internal_energy
    
    def set_value_keys(self):
        '''
        set the strings 
        '''
        self.pressure_key = "pressure"
        self.density_key = "density"
        self.velx_key = "vel_x"
        self.vely_key = "vel_y"
        self.vel_key = "vel_mag"
        self.inten_key = "int_energy"
        
        self.vel_str = "velocity"
        self.velx_str = "velocity (x)"
        self.vely_str = "velocity (y)"
        self.inten_str = "internal_energy"
        
    
    def run_check(self, ch_snap2):
        '''
        this is a temporary solution to checking whether two snaps come from the
            same ChollaRun instance. In the future, would like to incorporate
            class inheritance and have a instance called created_snaps for ChollaRun
            that will keep a list of created snaps
        
        for future, want to figure out how to verify that two children classes come
            from the same parent class
        
        temporary solution is to just ensure they have the same data directory
        '''
        err = False
        # err = (self.ch_snap.dataDir != ch_snap2.dataDir)
        
        # remove final snapshot number directory from dataDir string. compare both
        revdataDir = self.ch_snap.dataDir[::-1]
        lastslash_ind = revdataDir.index("/")
        datadir1 = revdataDir[lastslash_ind:][::-1]
        
        revdataDir = ch_snap2.dataDir[::-1]
        lastslash_ind = revdataDir.index("/")
        datadir2 = revdataDir[lastslash_ind:][::-1]
        err = (datadir1 != datadir2)
        
        if (err):
            err_message = f'''
            The two snapshots do not come from the same Cholla Run
            Currently can only compare snaps from same sim run
            '''
            
            raise Exception(err_message)
    
    def pressure(self, plt_kwargs=None):
        '''
        plots pressure, let users change plt kwargs for each plot if they want
        '''
        
        self.ch_snap.key_hydrodatacheck(self.pressure_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure(plt_kwargs)

    
    def density(self, plt_kwargs=None):
        '''
        plots density, let users change plt kwargs for each plot if they want
        '''

        self.ch_snap.key_hydrodatacheck(self.density_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.density(plt_kwargs)
    
    def internal_energy(self, plt_kwargs=None):
        '''
        plots internal energy, let users change plt kwargs for each plot if they want
        '''

        self.ch_snap.key_hydrodatacheck(self.inten_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.int_energy(plt_kwargs)
    
    def velocity_x(self, plt_kwargs=None):
        '''
        plots x velocity, let users change plt kwargs for each plot if they want
        '''
        
        
        self.ch_snap.key_hydrodatacheck(self.velx_key, key_str=self.velx_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_x(plt_kwargs)
    
    def velocity_y(self, plt_kwargs=None):
        '''
        plots y velocity, let users change plt kwargs for each plot if they want
        '''

        self.ch_snap.key_hydrodatacheck(self.vely_key, key_str=self.vely_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_y(plt_kwargs)
    
    def velocity(self, plt_kwargs=None):
        '''
        plots velocity, let users change plt kwargs for each plot if they want
        '''

        self.ch_snap.key_hydrodatacheck(self.vel_key, key_str=self.vel_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity(plt_kwargs)
        
        
    def density_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots density in comparison to ch_snap2's density
        '''

        
        self.ch_snap.key_hydrodatacheck(self.density_key, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.density_key, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.density_compare(ch_snap2, plt_kwargs)
        
    def pressure_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots pressure in comparison to ch_snap2's pressure
        '''

        self.ch_snap.key_hydrodatacheck(self.pressure_key, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.pressure_key, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure_compare(ch_snap2, plt_kwargs)
    
    def velocity_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots velocity in comparison to ch_snap2's velocity
        '''

        self.ch_snap.key_hydrodatacheck(self.vel_key, key_str=self.vel_str, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.vel_key, key_str=self.vel_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocity_compare(ch_snap2, plt_kwargs)
    
    def velocityx_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots x-velocity in comparison to ch_snap2's x-velocity
        '''
        
        self.ch_snap.key_hydrodatacheck(self.velx_key, key_str=self.velx_str, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.velx_key, key_str=self.velx_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityx_compare(ch_snap2, plt_kwargs)

    def velocityy_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots y-velocity in comparison to ch_snap2's y-velocity
        '''
        
        self.ch_snap.key_hydrodatacheck(self.vely_key, key_str=self.vely_str, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.vely_key, key_str=self.vely_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityy_compare(ch_snap2, plt_kwargs)
    
    def internalenergy_compare(self, ch_snap2, plt_kwargs=None):
        '''
        plots internal energy in comparison to ch_snap2's internal energy
        '''
        
        self.ch_snap.key_hydrodatacheck(self.inten_key, key_str=self.inten_str, raise_keyerr=True)
        ch_snap2.key_hydrodatacheck(self.inten_key, key_str=self.inten_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.internal_energy_compare(ch_snap2, plt_kwargs)

