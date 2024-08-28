import numpy as np
import h5py

class ChollaBoxHead:
    '''
    Cholla Box Head object
        Holds information regarding the location of the data

        Initialized with:
        - nBox (int): number of the box within snapshot
        - cell_offset_x (int): number of cells offset in x dimension 
        - cell_offset_y (int): number of cells offset in y dimension
        - cell_offset_z (int): number of cells offset in z dimension
        - local_dimx (int): number of cells in the box subvolume in x dimension
        - local_dimy (int): number of cells in the box subvolume in y dimension
        - local_dimz (int): number of cells in the box subvolume in z dimension
    '''
    
    def __init__(self, nBox, cell_offset_x, cell_offset_y, cell_offset_z, local_dimx, local_dimy, local_dimz):
        self.nBox = nBox
        self.offset = (cell_offset_x, cell_offset_y, cell_offset_z)
        self.local_dims = (local_dimx, local_dimy, local_dimz)
    
    def set_coords(self, xmin, ymin, zmin, dx, dy, dz):
        '''
        Set the domain coordinates based on some cell size and min
        
        Args:
            xmin (float): minimum x coordinate
            ymin (float): minimum y coordinate
            zmin (float): minimum z coordinate
            dx (float): cell size in x dimension
            dy (float): cell size in y dimension
            dz (float): cell size in z dimension
        Returns:
            ...
        '''

        # set the local min/maximum x coordinate of the box
        self.local_xmin = xmin + (self.offset[0] * dx)
        self.local_xmax = xmin + ((self.offset[0] + self.local_dims[0]) * dx)
        
        self.local_ymin = ymin + (self.offset[1] * dy)
        self.local_ymax = ymin + ((self.offset[1] + self.local_dims[1]) * dy)

        self.local_zmin = zmin + (self.offset[2] * dz)
        self.local_zmax = zmin + ((self.offset[2] + self.local_dims[2]) * dz)


        
class ChollaBox:
    '''
    Cholla Box object
        Holds information regarding the box for a specific snapshot 

        Initialized with:
            SnapPath (str): path to a snapshot directory
            chBoxHead (ChollaBoxHead): Cholla Box Head object connecting box to 
                rest of volume
            chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''
    
    def __init__(self, SnapPath, chBoxHead, chMacroFlags):
        self.BoxHead = ChollaBoxHead
        self.SnapPath = SnapPath

        # key strings hard coded in https://github.com/cholla-hydro/cholla/wiki/Output
        # using Makefile Parameter flags in https://github.com/cholla-hydro/cholla/wiki/Makefile-Parameters

        # hydro keys
        self.energy_str = "Energy"
        self.density_str = "density"
        self.momx_str = "momentum_x"
        self.momy_str = "momentum_y"
        self.momz_str = "momentum_z"
        self.gasenergy_str = "GasEnergy"
        
        if chMacroFlag.Cooling_GPU:
            # chemical species keys
            self.HIdensity_str = "HI_density"
            self.HIIdensity_str = "HII_density"
            self.HeIdensity_str = "HeI_density"
            self.HeIIdensity_str = "HeII_density"
            self.HeIIIdensity_str = "HeIII_density"
            self.edensity_str = "e_density"
            self.temp_str = "temperature"

            self.hydro_allkeys = {self.energy_str, self.density_str, self.momx_str, 
                                  self.momy_str, self.momz_str, self.gasenergy_str,
                                  self.HIdensity_str, self.HIIdensity_str,
                                  self.HeIdensity_str, self.HeIIdensity_str, self.HeIIIdensity_str,
                                  self.edensity_str, self.temp_str}
        else:
            self.hydro_allkeys = {self.energy_str, self.density_str, self.momx_str,
                                  self.momy_str, self.momz_str, self.gasenergy_str}


        if chMacroFlag.Particles:
            # particle keys
            self.densityCIC_str = "density"
            self.particleID_str = "particle_IDs"
            self.posx_str = "pos_x"
            self.posy_str = "pos_y"
            self.posz_str = "pos_z"
            self.velx_str = "vel_x"
            self.vely_str = "vel_y"
            self.velz_str = "vel_z"

            self.particle_allkeys = {self.densityCIC_str, self.particleID_str,
                                     self.posx_str, self.posy_str, self.posz_str,
                                     self.velx_str, self.vely_str, self.velz_str}
        else:
            self.particle_allkeys = {}
        
        if chMacroFlag.Gravity:
            # gravity keys
            self.potential_str = "potential"

            self.gravity_allkeys = {self.potential_str}
        else:
            self.gravity_allkeys = {}

    
    def get_nparts(self):
        '''
        Return the local number of particles in this box

        Args:
            ...
        Returns:
            ...
        '''

        fObj = h5py.File(self.get_particlefPath(), 'r')
        nparts = int(fObj.attrs['n_particles_local'])
        fObj.close()

        return nparts

    
    def get_hydrofPath(self):
        '''
        Return file path to hydro box file
        
        Args:
            ...
        Returns:
            fPath (str): file path to box's hydro data
        '''

        nSnap = self.SnapPath.split('/')[-1]
        fName = '{0}.{1}.{2}'.format(nSnap, 'h5', self.BoxHead.nBox)
        fPath = self.SnapPath + '/' + fName
        
        return fPath

    def get_particlefPath(self):
        '''
        Return file path to particle box file
        
        Args:
            ...
        Returns:
            fPath (str): file path to box's particle data
        '''

        nSnap = self.SnapPath.split('/')[-1]
        fName = '{0}_particles.{1}.{2}'.format(nSnap, 'h5', self.BoxHead.nBox)
        fPath = self.SnapPath + '/' + fName
        
        return fPath

    def get_gravityfPath(self):
        '''
        Return file path to gravity box file
        
        Args:
            ...
        Returns:
            fPath (str): file path to box's gravity data
        '''

        nSnap = self.SnapPath.split('/')[-1]
        fName = '{0}_gravity.{1}.{2}'.format(nSnap, 'h5', self.BoxHead.nBox)
        fPath = self.SnapPath + '/' + fName

        return fPath

    def check_hydrokey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed in hydro file

        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected hydro keys
        '''

        return data_key in self.hydro_allkeys

    def check_particlekey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed in particle file

        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected particle keys
        '''

        return data_key in self.particle_allkeys

    def check_gravitykey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed in gravity file

        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected gravity keys
        '''

        return data_key in self.gravity_allkeys

    def get_hydrodata(self, key, dtype=np.float32):
        '''
        Return a specific hydro dataset

        Args:
            key (str): key to access data from hydro hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        assert self.check_hydrokey(key)
        
        arr = np.zeros(self.BoxHead.local_dims, dtype=dtype)
        fObj = h5py.File(self.get_hydrofPath(), 'r')
        arr[:] = fObj.get(key)[:]
        fObj.close()

        return arr

    def get_particledata(self, key, dtype=np.float32):
        '''
        Return a specific particle dataset

        Args:
            key (str): key to access data from particle hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        assert self.check_particlekey(key)

        if (key == self.densityCIC_str):
            arr = np.zeros(self.BoxHead.local_dims, dtype=dtype)
        else:
            arr = np.zeros(self.get_nparts(), dtype=dtype)

        fObj = h5py.File(self.get_particlefPath(), 'r')
        arr[:] = fObj.get(key)[:]
        fObj.close()

        return arr

    def get_gravitydata(self, key, dtype=np.float32):
        '''
        Return a specific gravity dataset

        Args:
            key (str): key to access data from gravity hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        assert self.check_gravitykey(key)

        arr = np.zeros(self.BoxHead.local_dims, dtype=dtype)
        fObj = h5py.File(self.get_gravityfPath(), 'r')
        arr[:] = fObj.get(key)[:]
        fObj.close()

        return arr

    def place_data(self, data, arr):
        '''
        Place any data onto a larger array using box's offset
        
        Args:
            data (arr): data to be placed onto larger array
            arr (arr): will hold large data
        Returns:
            ...
        '''

        # TODO: assert that offset + local dims is less than arr size
        
        startX, startY, startZ = self.BoxHead.offset
        localX, localY, localZ = self.BoxHead.local_dims
        endX, endY, endZ = startX + localX, startY + localY, startZ + localZ

        arr[startX:endX, startY:endY, startZ:endZ] = data

