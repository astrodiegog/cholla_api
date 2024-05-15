import numpy as np
from time import time
from fast_histogram import histogram2d


class ChollaAnalysis:
    '''
    Cholla Analysis object
        Serves as the mediator between any large cool outputs requested and the
            calculations prescribed in ChollaCalculator. This object will use
            different outputs from calculator to complete a very specific task
            
    '''
    
    def __init__(self, dims, dtype=np.float32):
        self.dims = dims
        self.dtype = dtype
        
        # cgs constant values
        self.mp = 1.672622e-24
        self.kB = 1.380658e-16
        
        
    def create_arr(self):
        '''
        Create and return an empty array
        
        Args:
            ...
        Returns:
            (arr): array of initialized dimensions and datatype
        '''
        
        return np.zeros(self.dims, dtype=self.dtype)
    
    def create_subarr(self, new_dims):
        '''
        Create and return an empty array of new dimension shapes
        
        Args:
            new_dims (tuple): holds shape of new array
        Returns:
            (arr): array of initialized dimensions and datatype
        '''
        
        return np.zeros(new_dims, dtype=self.dtype)
    
    def create_hist2d(self, arr1, arr2, bins):
        '''
        Create a 2D histogram with well-defined bins.
        
        Args:
            arr1 (arr): array to be histogram-ed against arr2
            arr2 (arr): array to be histogram-ed against arr1
            bins (tuple): bins that can be passed onto numpy's histogram2d
                function (https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html)
        Returns:
            (tuple): 2D histogram, bin edges along 1st dimension, and bin edges 
                along 2nd dimension
        '''
        
        return np.histogram2d(arr1, arr2, bins=bins)
    
    def create_phase(self, log_temp, log_overdensity):
        '''
        Create a 2D histogram for temperature and overdensity. Designed for
            cosmological simulations
        
        Args:
            log_temp (arr): log-base10 of the gas temperature
            log_temp (arr): log-base10 of the gas overdensity
        Returns:
            (tuple): 2D histogram, bin edges along 1st dimension, and bin edges 
                along 2nd dimension
        '''
        min_loverdensity, max_loverdensity = -2, 4
        min_ltemp, max_ltemp = 2, 8
        
        log_overdensity_bins = np.linspace(min_loverdensity, max_loverdensity)
        log_temp_bins = np.linspace(min_ltemp, max_ltemp)
        
        # phase_bins = (log_overdensity_bins, log_temp_bins)
        # phasespace = self.create_hist2d(log_overdensity, log_temp, bins=phase_bins)
        
        phasespace = histogram2d(log_overdensity, log_temp, 
                                 range=((min_loverdensity, max_loverdensity), 
                                        (min_ltemp, max_ltemp)),
                                 bins=(49,49))
        
        return (phasespace, log_overdensity_bins, log_temp_bins)
    
    
    def stack_projection(self, k_index, k_domain, start, end, arr):
        '''
        Stack the density projection along a k_index 
            x-projection: k_index = 0
            y-projection: k_index = 1
            z-projection: k_index = 2
        
        Using starting and ending domain, place onto an array
        
        '''
        
        
        
    