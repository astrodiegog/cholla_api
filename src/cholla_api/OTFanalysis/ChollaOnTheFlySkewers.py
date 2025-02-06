import numpy as np
import h5py
from pathlib import Path

class ChollaOnTheFlySkewerHead:
    '''
    Cholla On The Fly Skewer Head

    Holds information regarding a specific individual skewer

        Initialized with:
        - skew_id (int): id of the skewer
        - n_i (int): length of the skewer
        - skew_key (str): string to access skewer

    '''
    def __init__(self, skew_id, n_i, skew_key):
        self.skew_id = skew_id
        self.n_i = n_i
        self.skew_key = skew_key


class ChollaOnTheFlySkewer:
    '''
    Cholla On The Fly Skewer
    
    Holds skewer specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFSkewerHead (ChollaOnTheFlySkewerHead): header
            information associated with skewer
        - fPath (PosixPath): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFSkewerHead, fPath):
        self.OTFSkewerHead = ChollaOTFSkewerHead
        self.fPath = fPath.resolve()    # convert to absolute path
        assert self.fPath.is_file()     # make sure file exists

        self.set_keys()                 # set possible skewer keys

    def set_keys(self):
        '''
        Check skewer group to set the available keys

        Args:
            ...
        Returns:
            ...
        '''

        keys_1D, keys_2D = [], []
        with h5py.File(self.fPath, 'r') as fObj:
            self.allkeys = set(fObj[self.OTFSkewersiHead.skew_key].keys())
            for key in self.allkeys:
                if fObj[self.OTFSkewersiHead.skew_key].get(key).ndim == 1:
                    keys_1D.append(key)
                if fObj[self.OTFSkewersiHead.skew_key].get(key).ndim == 2:
                    keys_2D.append(key)

        self.keys_1D = set(keys_1D)
        self.keys_2D = set(keys_2D)

        return

    def check_datakey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed in skewers file

        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected data keys
        '''

        return data_key in self.allkeys

    def get_skewerdata(self, key, dtype=np.float32):
        '''
        Return a specific skewer dataset

        Args:
            key (str): key to access data from hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        assert self.check_datakey(key)

        arr = np.zeros(self.OTFSkewerHead.n_i, dtype=dtype)
        with h5py.File(self.fPath, 'r') as fObj:
            arr[:] = fObj[self.OTFSkewerHead.skew_key].get(key)[self.OTFSkewerHead.skew_id, :]

        return arr

    def get_density(self, dtype=np.float32):
        '''
        Return the density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): density
        '''
        density_str = 'density'

        return self.get_skewerdata(density_str, dtype=dtype)

    def get_HIdensity(self, dtype=np.float32):
        '''
        Return the HI density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HI density
        '''
        HI_str = 'HI_density'

        return self.get_skewerdata(HI_str, dtype=dtype)

    def get_HeIIdensity(self, dtype=np.float32):
        '''
        Return the HeII density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HeII density
        '''
        HeII_str = 'HeII_density'

        return self.get_skewerdata(HeII_str, dtype=dtype)

    def get_losvelocity(self, dtype=np.float32):
        '''
        Return the line-of-sight velocity array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): line-of-sight velocity
        '''
        vel_str = 'los_velocity'

        return self.get_skewerdata(vel_str, dtype=dtype)

    def get_temperature(self, dtype=np.float32):
        '''
        Return the temperature array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): temperature
        '''
        temp_str = 'temperature'

        return self.get_skewerdata(temp_str, dtype=dtype)




class ChollaOnTheFlySkewers_iHead:
    '''
    Cholla On The Fly Skewers_i Head

    Holds information regarding a specific skewer hdf5 group

        Initialized with:
        - n_i (int): length of the skewers
        - n_j (int): length of first dimension spanning cube
        - n_k (int): lenth of second dimension spanning cube
        - n_stride (int): stride cell number between skewers
        - skew_key (str): string to access skewer
    '''
    def __init__(self, n_i, n_j, n_k, n_stride, skew_key):
        self.n_i = n_i
        self.n_j = n_j
        self.n_k = n_k
        self.n_stride = n_stride
        self.skew_key = skew_key

        # number of skewers, assumes nstride is same along both j and k dims
        self.n_skews = int( (self.n_j * self.n_k) / (self.n_stride * self.n_stride) )


class ChollaOnTheFlySkewers_i:
    '''
    Cholla On The Fly Skewers
    
    Holds skewer specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFSkewersiHead (ChollaOnTheFlySkewers_iHead): header
            information associated with skewer
        - fPath (PosixPath): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFSkewersiHead, fPath):
        self.OTFSkewersiHead = ChollaOTFSkewersiHead
        self.fPath = fPath.resolve()    # convert to absolute path
        assert self.fPath.is_file()     # make sure file exists

        self.set_keys()                 # set possible skewer keys

    def set_keys(self):
        '''
        Check skewer group to set the available keys

        Args:
            ...
        Returns:
            ...
        '''

        keys_1D, keys_2D = [], []
        with h5py.File(self.fPath, 'r') as fObj:
            self.allkeys = set(fObj[self.OTFSkewersiHead.skew_key].keys())
            for key in self.allkeys:
                if fObj[self.OTFSkewersiHead.skew_key].get(key).ndim == 1:
                    keys_1D.append(key)
                if fObj[self.OTFSkewersiHead.skew_key].get(key).ndim == 2:
                    keys_2D.append(key)

        self.keys_1D = set(keys_1D)
        self.keys_2D = set(keys_2D)

        return

    def check_datakey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed in skewers file

        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected data keys
        '''

        return data_key in self.allkeys

    def get_skeweralldata(self, key, dtype=np.float32):
        '''
        Return a specific dataset for all skewers.
            Use this method with caution, as the resulting array can be large

            For (2048)^3 + nstride=4 + float64, resulting array will be ~4 GBs

        Args:
            key (str): key to access data from hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        assert self.check_datakey(key)

        if key in self.keys_1D:
            arr = np.zeros((self.OTFSkewersiHead.n_skews), dtype=dtype)
            with h5py.File(self.fPath, 'r') as fObj:
                arr[:] = fObj[self.OTFSkewersiHead.skew_key].get(key)[:]
        elif key in self.keys_2D:
            arr = np.zeros((self.OTFSkewersiHead.n_skews, self.OTFSkewersiHead.n_i), dtype=dtype)
            with h5py.File(self.fPath, 'r') as fObj:
                arr[:,:] = fObj[self.OTFSkewersiHead.skew_key].get(key)[:, :]

        return arr

    def get_alldensity(self, dtype=np.float32):
        '''
        Return density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): density
        '''
        density_str = 'density'

        return self.get_skeweralldata(density_str, dtype=dtype)

    def get_allHIdensity(self, dtype=np.float32):
        '''
        Return HI density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HI density
        '''
        HI_str = 'HI_density'

        return self.get_skeweralldata(HI_str, dtype=dtype)

    def get_allHeIIdensity(self, dtype=np.float32):
        '''
        Return HeII density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): phase space
        '''
        HeII_str = 'HeII_density'

        return self.get_skeweralldata(HeII_str, dtype=dtype)

    def get_alllosvelocity(self, dtype=np.float32):
        '''
        Return line-of-sight velocity array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): line-of-sight velocity
        '''
        vel_str = 'los_velocity'

        return self.get_skeweralldata(vel_str, dtype=dtype)

    def get_alltemperature(self, dtype=np.float32):
        '''
        Return temperature array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): temperature
        '''
        temp_str = 'temperature'

        return self.get_skeweralldata(temp_str, dtype=dtype)

    def get_skewer_obj(self, skewid):
        '''
        Return ChollaOnTheFlySkewer object of this analysis

        Args:
            skew_id (int): skewer id
        Return:
            OTFSkewer (ChollaOnTheFlySkewer): skewer object
        '''
        OTFSkewerHead = ChollaOnTheFlySkewerHead(skewid, self.OTFSkewersiHead.n_i, 
                                                 self.OTFSkewersiHead.skew_key)

        return ChollaOnTheFlySkewer(OTFSkewerHead, self.fPath)

    def get_dvHubble(self, dvHubble):
        '''
        Return the Hubble flow velocity along line-of-sight of skewers

        Args:
            dvHubble (float): differential Hubble flow velocity
        Returns:
            (arr): Hubble velocity along skewer line-of-sight
        '''

        iter_arr = np.arange(self.OTFSkewersiHead.n_i)

        return dvHubble * (iter_arr + 0.5)




class ChollaOnTheFlySkewers:
    '''
    Cholla On The Fly Skewers
    
    Holds on-the-fly skewers specific information to an output with methods to 
            access data for that output

        Initialized with:
        - fPath (PosixPath): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, nSkewer, SkewersPath):
        self.OTFSkewersfPath = fPath.resolve() # convert to absolute path
        assert self.OTFSkewersfPath.is_file() # make sure file exists

        self.xskew_str = "skewers_x"
        self.yskew_str = "skewers_y"
        self.zskew_str = "skewers_z"

        # set grid information (ncells, dist between cells, nstride)
        self.set_gridinfo()
        dx_Mpc = self.dx / 1.e3 # [Mpc]
        dy_Mpc = self.dy / 1.e3
        dz_Mpc = self.dz / 1.e3

        # set cosmology params
        self.set_cosmoinfo()

        # grab current hubble param & info needed to calculate hubble flow
        H = self.get_currH()
        cosmoh = self.H0 / 100.

        # calculate proper distance along each direction
        dxproper = dx_Mpc * self.current_a / cosmoh # [h-1 Mpc]
        dyproper = dy_Mpc * self.current_a / cosmoh
        dzproper = dz_Mpc * self.current_a / cosmoh

        # calculate hubble flow through a cell along each axis
        self.dvHubble_x = H * dxproper # [km s-1]
        self.dvHubble_y = H * dyproper
        self.dvHubble_z = H * dzproper

    def set_gridinfo(self, datalength_str='density'):
        '''
        Set grid information by looking at attribute of file object and shape of 
            data sets
        
        Args:
            - datalength_str (str): (optional) key to dataset used to find the
                number of skewers and cells along an axis
        Returns:
            ...
        '''

        with h5py.File(self.OTFSkewersfPath, 'r') as fObj:
            # grab length of box in units of [kpc]
            Lx, Ly, Lz = np.array(fObj.attrs['Lbox'])

            # set number of skewers and stride number along each direction 
            self.nskewersx, self.nx = fObj[self.xskew_str][datalength_str].shape
            self.nskewersy, self.ny = fObj[self.yskew_str][datalength_str].shape
            self.nskewersz, self.nz = fObj[self.zskew_str][datalength_str].shape


        # we know nskewers_i = (nj * nk) / (nstride_i * nstride_i)
        # so nstride_i = sqrt( (nj * nk) / (nskewers_i) )
        self.nstride_x = int(np.sqrt( (self.ny * self.nz)/(self.nskewersx) ))
        self.nstride_y = int(np.sqrt( (self.nz * self.nx)/(self.nskewersy) ))
        self.nstride_z = int(np.sqrt( (self.nx * self.ny)/(self.nskewersz) ))

        # save cell distance in each direction to later calculate hubble flow
        self.dx = Lx / self.nx # [kpc]
        self.dy = Ly / self.ny
        self.dz = Lz / self.nz

    def set_cosmoinfo(self):
        '''
        Set cosmological attributes for this object

        Args:
            ChollaOTFAnalysis (ChollaOnTheFlyAnalysis): specific Analysis object to use
        Returns:
            ...
        '''

        with h5py.File(self.OTFSkewersfPath, 'r') as fObj:
            self.Omega_R = fObj.attrs['Omega_R'].item()
            self.Omega_M = fObj.attrs['Omega_M'].item()
            self.Omega_L = fObj.attrs['Omega_L'].item()
            self.Omega_K = fObj.attrs['Omega_K'].item()

            self.w0 = fObj.attrs['w0'].item()
            self.wa = fObj.attrs['wa'].item()

            self.H0 = fObj.attrs['H0'].item() # expected in km/s/Mpc
            self.current_a = fObj.attrs['current_a'].item()
            self.current_z = fObj.attrs['current_z'].item()

        return

    def get_currH(self):
        '''
        Return the Hubble parameter at the current scale factor

        Args:
            ...
        Returns:
            H (float): Hubble parameter (km/s/Mpc)
        '''

        a2 = self.current_a * self.current_a
        a3 = a2 * self.current_a
        a4 = a3 * self.current_a
        DE_factor = (self.current_a)**(-3. * (1. + self.w0 + self.wa))
        DE_factor *= np.exp(-3. * self.wa * (1. - self.current_a))

        H0_factor = (self.Omega_R / a4) + (self.Omega_M / a3)
        H0_factor += (self.Omega_K / a2) + (self.Omega_L * DE_factor)

        return self.H0 * np.sqrt(H0_factor)

    def get_skewersx_obj(self):
        '''
        Return ChollaOnTheFlySkewers_i object of the x-skewers

        Args:
            ...
        Return:
            OTFSkewerx (ChollaOnTheFlySkewers_i): skewer object
        '''

        OTFSkewersxHead = ChollaOnTheFlySkewers_iHead(self.nx, self.ny, self.nz,
                                                      self.nstride_x, self.xskew_str)

        OTFSkewerx = ChollaOnTheFlySkewers_i(OTFSkewersxHead, self.OTFSkewersfPath)

        return OTFSkewerx

    def get_skewersy_obj(self):
        '''
        Return ChollaOnTheFlySkewers_i object of the y-skewers

        Args:
            ...
        Return:
            OTFSkewery (ChollaOnTheFlySkewers_i): skewer object
        '''

        OTFSkewersyHead = ChollaOnTheFlySkewers_iHead(self.ny, self.nx, self.nz,
                                                      self.nstride_y, self.yskew_str)

        OTFSkewery = ChollaOnTheFlySkewers_i(OTFSkewersyHead, self.OTFSkewersfPath)

        return OTFSkewery

    def get_skewersz_obj(self):
        '''
        Return ChollaOnTheFlySkewers_i object of the z-skewers

        Args:
            ...
        Return:
            OTFSkewerz (ChollaOnTheFlySkewers_i): skewer object
        '''

        OTFSkewerszHead = ChollaOnTheFlySkewers_iHead(self.nz, self.nx, self.ny,
                                                      self.nstride_z, self.zskew_str)

        OTFSkewerz = ChollaOnTheFlySkewers_i(OTFSkewerszHead, self.OTFSkewersfPath)

        return OTFSkewerz




    
