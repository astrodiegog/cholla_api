import numpy as np


class ChollaExpansionHistory:
    '''
    Cholla Expansion History object
        Holds the history of a cosmology run
            
        Initialized with:
            expansionfPath (str): path to a directory holding all snapshot directories
    '''

    def __init__(self, expansion_fname):
        self.fPath = expansion_fname

    def get_time(self):
        '''
        Return the time of an expansion history

        Args:
            ...
        Returns:
            arr (arr): time in Myrs
        '''
        arr = np.loadtxt(self.fPath)

        return arr[:,0]

    def get_scalefactor(self):
        '''
        Return the scale factor of an expansion history

        Args:
            ...
        Returns:
            arr (arr): scale factor
        '''
        arr = np.loadtxt(self.fPath)

        return arr[:,1]



