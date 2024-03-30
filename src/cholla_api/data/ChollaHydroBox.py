import numpy as np
import h5py

class ChollaHydroBoxHead:
    '''
    Cholla Hydro Data Box Head object
        Holds information regarding the hydro data for a specific snapshot 
            within a specific box
        Initialized with:
        - nBox (int): number of the box within snapshot
    '''
    
    def __init__(self, nBox):
        self.nBox = nBox
        self.head_set = False
        
    def set_head(self, nSnap, namebase, dataDir, old_format):
        '''
        Set the header attributes for this object.
        
        Args:
            nSnap (int): what snapshot number to use
            nBox (int): what box number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        if nSnap == 0:
            return
        fName = '{0}.{1}.{2}'.format(nSnap, namebase, self.nBox)
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(nSnap) + '/' + fName
        fObj = h5py.File(fPath, 'r')
        
        self.offset = fObj.attrs['offset']
        self.local_dims = fObj.attrs['dims_local']
        
        fObj.close()
        
        self.head_set = True

        
class ChollaHydroBox:
    '''
    Cholla Hydro Box object
        Holds information regarding the hydro data for a specific snapshot 
            within a specific box.
        Contains methods to 
            get box's file path 
            get a specific data type
            place a specific data type onto an array
    '''
    
    def __init__(self, HydroBoxHead, nSnap):
        self.head = HydroBoxHead
        assert self.head.head_set
        self.nBox = self.head.nBox
        self.nSnap = nSnap
    
    def get_fPath(self, namebase, dataDir, old_format):
        '''
        Return file path to this box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            fPath (str): file path to this box
        '''
        fName = '{0}.{1}.{2}'.format(self.nSnap, namebase, self.nBox)
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(self.nSnap) + '/' + fName
        return fPath
    
    def get_data(self, namebase, dataDir, data_key, old_format):
        '''
        Return a specific dataset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            old_format (bool): whether the old file structure was used
        Returns:
            arr (arr): dataset
        '''
        
        arr = np.zeros(self.head.local_dims)
        fObj = h5py.File(self.get_fPath(namebase, dataDir, old_format))
        arr[:] = np.array(fObj[data_key]).reshape(self.head.local_dims)
        fObj.close()
        
        return arr
    
    def place_data(self, namebase, dataDir, data_key, arr, old_format):
        '''
        Place a specific dataset onto a larger array using box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            arr (arr): array to place dataset
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        startX, startY, startZ = self.head.offset
        localX, localY, localZ = self.head.local_dims
        endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

        arr[startX:endX, startY:endY, startZ:endZ] = self.get_data(namebase, 
                                                                   dataDir, 
                                                                   data_key,
                                                                   old_format)
        
        
    def place_anydata(self, data, arr):
        '''
        Place any data onto a larger array using box's offset
        
        Args:
            data (arr): data to be placed onto larger array
            arr (arr): will hold large data
        Returns:
            ...
        '''
        
        startX, startY, startZ = self.head.offset
        localX, localY, localZ = self.head.local_dims
        endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

        arr[startX:endX, startY:endY, startZ:endZ] = data

