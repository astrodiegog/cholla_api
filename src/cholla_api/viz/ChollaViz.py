'''
Handles all visualization API with user. Should be simple
'''


'''
when catching the keyerror for each value, should I be raising
    KeyError or just returning a None?
'''

from cholla_api.viz.Cholla2DViz import Cholla2DViz

class ChollaViz:
    
    def __init__(self, chollasnap, num_D, test_name="", plt_kwargs={}):
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
    
    def pressure(self, plt_kwargs=None):
        '''
        let users change plt kwargs for each plot if they want
        '''
        
        # make sure we have pressure before plotting, if not, raise keyerr
        try:
            self.ch_snap.data["pressure"]
        except KeyError:
            err_message = '-- Error --\n'
            err_message += 'Unable to plot pressure \n'
            err_message += 'Missing pressure in snapshot'
            print(err_message)
            raise
        
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
        
        self.chollaxDViz.pressure(plt_kwargs)
    
    def density(self, plt_kwargs=None):
        
        # make sure we have density before plotting, if not, raise keyerr
        try:
            self.ch_snap.data["density"]
        except KeyError:
            err_message = '-- Error --'
            err_message = 'Unable to plot density \n'
            err_message += 'Missing density in snapshot'
            print(err_message)
            raise
            
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.density(plt_kwargs)
    
    def velocity_x(self, plt_kwargs=None):
        
        # make sure we have velocity_x before plotting, if not, raise keyerr
        try:
            self.ch_snap.data["vel_x"]
        except KeyError:
            err_message = '-- Error --'
            err_message = 'Unable to plot velocity_x \n'
            err_message += 'Missing vel_x in snapshot'
            print(err_message)
            raise
            
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_x(plt_kwargs)
    
    def velocity_y(self, plt_kwargs=None):
        
        # make sure we have velocity_y before plotting, if not, raise keyerr
        try:
            self.ch_snap.data["vel_y"]
        except KeyError:
            err_message = '-- Error --'
            err_message = 'Unable to plot velocity_y \n'
            err_message += 'Missing vel_y in snapshot'
            print(err_message)
            raise
            
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity_y(plt_kwargs)
    
    def velocity(self, plt_kwargs=None):
        
        # make sure we have velocity before plotting, if not, raise keyerr
        try:
            self.ch_snap.data["vel_mag"]
        except KeyError:
            err_message = '-- Error --'
            err_message = 'Unable to plot velocity \n'
            err_message += 'Missing vel_mag in snapshot'
            print(err_message)
            raise
            
        if plt_kwargs is None:
            plt_kwargs = self.plt_kwargs
            
        self.chollaxDViz.velocity(plt_kwargs)
        
        
