'''
Handles visualization comparison between snapshots.
    this is essentially a wrapper for chollaviz.value_compare()
'''



import numpy as np

from cholla_api.viz.ChollaxDViz import ChollaxDViz


class ChollaVizCompare:
    '''
    Class that compares two snapshots, very similar to ChollaViz, but just for comparing snapshots
    '''
    def __init__(self, ch_snap1, ch_snap2, test_name="", plt_kwargs={"show": True, "save": False}, reduce_dim_fn=None):
        self.ch_snap1 = ch_snap1
        self.ch_snap2 = ch_snap2
        self.test_name = test_name
        self.plt_kwargs = plt_kwargs
        self.reduce_dim_fn = reduce_dim_fn
        self.num_D = int(np.sum(self.ch_snap1.head['dims'] != 1.))
        
        self.check_dims(raise_err=True)
        
        if self.num_D == 1:
            self.chollaxDViz = ChollaxDViz(self.ch_snap1, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
            
        if self.num_D == 2:
            self.chollaxDViz = ChollaxDViz(self.ch_snap1, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
            
        if self.num_D == 3:
            self.chollaxDViz = ChollaxDViz(self.ch_snap1, num_dims=self.num_D, test_name=self.test_name, reduce_dim_fn=self.reduce_dim_fn)
    

    def check_dims(self, raise_err=False):
        '''
        make sure both snaps have the same dimensions
        '''

        try:
            assert np.array_equal(self.ch_snap1.head["dims"], self.ch_snap2.head["dims"])
        except AssertionError:
            err_message = '-- Error --\n'
            err_message += f'Snapshots have different dimensions, unable to compare'
            print(err_message)
            if raise_keyerr:
                raise    


    
    def density(self, plt_kwargs=None):
        
        density_key = "density"
        
        self.ch_snap1.key_datacheck(density_key, raise_keyerr=True)
        self.ch_snap2.key_datacheck(density_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.density_compare(self.ch_snap2, plt_kwargs)
        
    def pressure(self, plt_kwargs=None):
        
        pressure_key = "pressure"
        
        self.ch_snap1.key_datacheck(pressure_key, raise_keyerr=True)
        self.ch_snap2.key_datacheck(pressure_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure_compare(self.ch_snap2, plt_kwargs)
    
    def velocity(self, plt_kwargs=None):
        
        vel_key = "vel_mag"
        vel_str = "velocity"
        
        self.ch_snap1.key_datacheck(vel_key, key_str=vel_str, raise_keyerr=True)
        self.ch_snap2.key_datacheck(vel_key, key_str=vel_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocity_compare(self.ch_snap2, plt_kwargs)
    
    def velocityx(self, plt_kwargs=None):
        
        velx_key = "vel_x"
        velx_str = "velocity (x)"
        
        self.ch_snap1.key_datacheck(velx_key, key_str=velx_str, raise_keyerr=True)
        self.ch_snap2.key_datacheck(velx_key, key_str=velx_str, raise_keyerr=True)
        
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityx_compare(self.ch_snap2, plt_kwargs)

    def velocityy(self, ch_snap2, plt_kwargs=None):
        
        vely_key = "vel_y"
        vely_str = "velocity (y)"
        
        self.ch_snap1.key_datacheck(vely_key, key_str=vely_str, raise_keyerr=True)
        self.ch_snap2.key_datacheck(vely_key, key_str=vely_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityy_compare(self.ch_snap2, plt_kwargs)

