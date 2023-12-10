'''
Handles all visualization API with user. Should be simple
'''


'''
when catching the keyerror for each value, should I be raising
    KeyError or just returning a None?
'''

from cholla_api.viz.Cholla1DViz import Cholla1DViz
from cholla_api.viz.Cholla2DViz import Cholla2DViz

class ChollaViz:
    
    def __init__(self, chollasnap, num_D, test_name="", plt_kwargs={"show": True, "save": False}):
        self.ch_snap = chollasnap
        self.test_name = test_name
        self.num_D = num_D
        self.plt_kwargs = plt_kwargs
        
        if self.num_D == 1:
            self.chollaxDViz = Cholla1DViz(self.ch_snap, test_name)
        if self.num_D == 2:
            self.chollaxDViz = Cholla2DViz(self.ch_snap, test_name)
        if self.num_D == 3:
            self.chollaxDViz = Cholla3DViz(self.ch_snap, test_name)
    
    
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
        
        
        if (self.ch_snap.dataDir != ch_snap2.dataDir):
            err_message = f'''
            The two snapshots do not come from the same Cholla Run
            Currently can only compare snaps from same sim run
            '''
            
            raise Exception(err_message)
    
    def pressure(self, plt_kwargs=None):
        '''
        let users change plt kwargs for each plot if they want
        '''
        
        pressure_key = "pressure"
        self.ch_snap.key_datacheck(pressure_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure(plt_kwargs)
    
    def density(self, plt_kwargs=None):
            
        density_key = "density"
        self.ch_snap.key_datacheck(density_key, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.density(plt_kwargs)
    
    def velocity_x(self, plt_kwargs=None):
        
        velx_key = "vel_x"
        velx_str = "velocity (x)"
        self.ch_snap.key_datacheck(velx_key, key_str=velx_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_x(plt_kwargs)
    
    def velocity_y(self, plt_kwargs=None):
        
        vely_key = "vel_y"
        vely_str = "velocity (y)"
        self.ch_snap.key_datacheck(vely_key, key_str=vely_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_y(plt_kwargs)
    
    def velocity(self, plt_kwargs=None):
        
        vel_key = "vel_mag"
        vel_str = "velocity"
        self.ch_snap.key_datacheck(vel_key, key_str=vel_str, raise_keyerr=True)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity(plt_kwargs)
        
        
    def density_compare(self, ch_snap2, plt_kwargs=None):
        
        density_key = "density"
        
        self.ch_snap.key_datacheck(density_key, raise_keyerr=True)
        ch_snap2.key_datacheck(density_key, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.density_compare(ch_snap2, plt_kwargs)
        
    def pressure_compare(self, ch_snap2, plt_kwargs=None):
        
        pressure_key = "pressure"
        
        self.ch_snap.key_datacheck(pressure_key, raise_keyerr=True)
        ch_snap2.key_datacheck(pressure_key, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure_compare(ch_snap2, plt_kwargs)
    
    def velocity_compare(self, ch_snap2, plt_kwargs=None):
        
        vel_key = "vel_mag"
        vel_str = "velocity"
        
        self.ch_snap.key_datacheck(vel_key, key_str=vel_str, raise_keyerr=True)
        ch_snap2.key_datacheck(vel_key, key_str=vel_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocity_compare(ch_snap2, plt_kwargs)
    
    def velocityx_compare(self, ch_snap2, plt_kwargs=None):
        
        velx_key = "vel_x"
        velx_str = "velocity (x)"
        
        self.ch_snap.key_datacheck(velx_key, key_str=velx_str, raise_keyerr=True)
        ch_snap2.key_datacheck(velx_key, key_str=velx_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityx_compare(ch_snap2, plt_kwargs)

    def velocityy_compare(self, ch_snap2, plt_kwargs=None):
        
        vely_key = "vel_y"
        vely_str = "velocity (y)"
        
        self.ch_snap.key_datacheck(vely_key, key_str=vely_str, raise_keyerr=True)
        ch_snap2.key_datacheck(vely_key, key_str=vely_str, raise_keyerr=True)
        
        self.run_check(ch_snap2)
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.velocityy_compare(ch_snap2, plt_kwargs)