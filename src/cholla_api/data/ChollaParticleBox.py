import numpy as np
import h5py

class ChollaParticleBoxHead:
    '''
    Cholla Particle Data Box Head object
        Holds information regarding the particle data for a specific snapshot 
            within a specific box
        Initialized with:
        - nBox (int): number of the box within snapshot
    '''
    
    def __init__(self, nBox):
        self.nBox = nBox
        self.head_set = False
    
    def set_head(self, nSnap, nBox, namebase, dataDir, old_format):
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
        fName = '{0}_particles.{1}.{2}'.format(nSnap, namebase, nBox)
        
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(nSnap) + '/' + fName
        
        fObj = h5py.File(fPath, 'r')
        
        self.offset = fObj.attrs['offset']
        self.local_dims = fObj.attrs['dims_local']
        
        self.local_nparts = int(fObj.attrs['n_particles_local'])
        
        fObj.close()
        
        self.head_set = True


        
class ChollaParticleBox:
    '''
    Cholla Particle Box object
        Holds information regarding the particle data for a specific snapshot 
            within a specific box.
        Contains methods to 
            get box's file path 
            get a specific data type
            place a specific data type onto an array
    '''
    
    def __init__(self, ParticleBoxHead, nSnap):
        self.head = ParticleBoxHead
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
        fName = '{0}_particles.{1}.{2}'.format(self.nSnap, namebase, self.nBox)
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
        density_key = "density" # density cloud-in-cell
        
        fObj = h5py.File(self.get_fPath(namebase, dataDir, old_format))
        if (data_key == density_key):
            arr = np.zeros(self.head.local_dims)
            arr[:] = np.array(fObj[data_key]).reshape(self.head.local_dims)
        else:
            arr = np.zeros(self.head.local_nparts)
            arr[:] = np.array(fObj[data_key]).reshape(self.head.local_nparts)
        fObj.close()
        
        return arr
    
  
    def place_data(self, namebase, dataDir, data_key, arr, old_format):
        '''
        Place a specific dataset onto a larger array
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            arr (arr): array to place dataset
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        density_key = "density" # density cloud-in-cell
        if (data_key == density_key):
            startX, startY, startZ = self.head.offset
            localX, localY, localZ = self.head.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            arr[startX:endX, startY:endY, startZ:endZ] = self.get_data(namebase, 
                                                                       dataDir, 
                                                                       data_key, 
                                                                       old_format)
        else:
            startN = self.head.offset_nparts
            localN = self.head.local_nparts
            endN = startN + localN
            
            arr[startN:endN] = self.get_data(namebase, dataDir, data_key, 
                                             old_format)
