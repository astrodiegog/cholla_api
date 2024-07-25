import numpy as np
import h5py


class ChollaOnTheFlySkewerHead:
    '''
    Cholla On The Fly Skewer Head

    Holds information regarding a specific individual skewer

        Initialized with:
        - skew_id (int): id of the skewer
        - n_i (int): length of the skewer
        - skew_key (str): string to access skewer

    Values are returned in code units unless otherwise specified.
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
        - fPath (str): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFSkewerHead, fPath):
        self.OTFSkewerHead = ChollaOTFSkewerHead
        self.fPath = fPath

        self.HI_str = 'HI_density'
        self.HeII_str = 'HeII_density'
        self.density_str = 'density'
        self.vel_str = 'los_velocity'
        self.temp_str = 'temperature'

        self.allkeys = {self.HI_str, self.HeII_str, self.density_str,
                        self.vel_str, self.temp_str}

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
        fObj = h5py.File(self.fPath, 'r')
        arr[:] = fObj[self.OTFSkewerHead.skew_key].get(key)[self.OTFSkewerHead.skew_id, :]
        fObj.close()

        return arr

    def get_density(self, dtype=np.float32):
        '''
        Return the density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): density
        '''

        return self.get_skewerdata(self.density_str, dtype=dtype)

    def get_HIdensity(self, dtype=np.float32):
        '''
        Return the HI density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HI density
        '''

        return self.get_skewerdata(self.HI_str, dtype=dtype)

    def get_HeIIdensity(self, dtype=np.float32):
        '''
        Return the HeII density array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HeII density
        '''

        return self.get_skewerdata(self.HeII_str, dtype=dtype)

    def get_losvelocity(self, dtype=np.float32):
        '''
        Return the line-of-sight velocity array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): line-of-sight velocity
        '''

        return self.get_skewerdata(self.vel_str, dtype=dtype)

    def get_temperature(self, dtype=np.float32):
        '''
        Return the temperature array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): temperature
        '''

        return self.get_skewerdata(self.temp_str, dtype=dtype)




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

    Values are returned in code units unless otherwise specified.
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
        - fPath (str): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFSkewersiHead, fPath):
        self.OTFSkewersiHead = ChollaOTFSkewersiHead
        self.fPath = fPath

        self.HI_str = 'HI_density'
        self.HeII_str = 'HeII_density'
        self.density_str = 'density'
        self.vel_str = 'los_velocity'
        self.temp_str = 'temperature'

        self.allkeys = {self.HI_str, self.HeII_str, self.density_str,
                        self.vel_str, self.temp_str}

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

        arr = np.zeros((self.OTFSkewersiHead.n_skews, self.OTFSkewersiHead.n_i), dtype=dtype)
        fObj = h5py.File(self.fPath, 'r')
        arr[:,:] = fObj[self.OTFSkewersiHead.skew_key].get(key)[:, :]
        fObj.close()

        return arr

    def get_alldensity(self, dtype=np.float32):
        '''
        Return density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): density
        '''

        return self.get_skeweralldata(self.density_str, dtype=dtype)

    def get_allHIdensity(self, dtype=np.float32):
        '''
        Return HI density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): HI density
        '''

        return self.get_skeweralldata(self.HI_str, dtype=dtype)

    def get_allHeIIdensity(self, dtype=np.float32):
        '''
        Return HeII density array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): phase space
        '''

        return self.get_skeweralldata(self.HeII_str, dtype=dtype)

    def get_alllosvelocity(self, dtype=np.float32):
        '''
        Return line-of-sight velocity array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): line-of-sight velocity
        '''

        return self.get_skeweralldata(self.vel_str, dtype=dtype)

    def get_alltemperature(self, dtype=np.float32):
        '''
        Return temperature array for all skewers

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): temperature
        '''

        return self.get_skeweralldata(self.temp_str, dtype=dtype)

    def get_skewer_obj(self, skewid):
        '''
        Return ChollaOnTheFlySkewer object of this analysis

        Args:
            skew_id (int): skewer id
        Return:
            OTFSkewer (ChollaOnTheFlySkewer): skewer object
        '''
        OTFSkewerHead = ChollaOnTheFlySkewerHead(skew_id, self.OTFSkewersiHead.n_i, 
                                                 self.ChollaOTFSkewersiHead.skew_key)

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
        - nSkewer (nSkewer): number of the skewer output
        - SkewersPath (str): directory path to skewer output files
        - ChollaGrid (ChollaGrid): grid holding domain information

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, nSkewer, SkewersPath, ChollaGrid):
        self.OTFSkewersfPath = SkewersPath + '/' + str(nSkewer) + '_skewers.h5'

        # current implementation only works for cube grid
        assert ChollaGrid.nx_global == ChollaGrid.ny_global
        assert ChollaGrid.nx_global == ChollaGrid.nz_global
        assert ChollaGrid.Lx == ChollaGrid.Ly
        assert ChollaGrid.Lx == ChollaGrid.Lz

        self.nx = ChollaGrid.nx_global
        self.ny = ChollaGrid.ny_global
        self.nz = ChollaGrid.nz_global

        dx = ChollaGrid.dx / 1e3 # convert kpc --> Mpc

        # set cosmology params
        self.set_cosmoinfo()

        # grab current hubble param & info needed to calculate hubble flow
        H = self.get_currH()
        cosmoh = self.H0 / 100.
        dxproper = dx * self.current_a / cosmoh

        # calculate hubble flow through a cell
        self.dvHubble = H * dxproper

        self.xskew_str = "skewers_x"
        self.yskew_str = "skewers_y"
        self.zskew_str = "skewers_z"

    def set_cosmoinfo(self):
        '''
        Set cosmological attributes for this object

        Args:
            ChollaOTFAnalysis (ChollaOnTheFlyAnalysis): specific Analysis object to use
        Returns:
            ...
        '''

        fObj = h5py.File(self.OTFSkewersfPath, 'r')

        self.Omega_R = fObj.attrs['Omega_R'].item()
        self.Omega_M = fObj.attrs['Omega_M'].item()
        self.Omega_L = fObj.attrs['Omega_L'].item()
        self.Omega_K = fObj.attrs['Omega_K'].item()

        self.w0 = fObj.attrs['w0'].item()
        self.wa = fObj.attrs['wa'].item()

        self.H0 = fObj.attrs['H0'].item() # expected in km/s/Mpc
        self.current_a = fObj.attrs['current_a'].item()
        self.current_z = fObj.attrs['current_z'].item()

        fObj.close()

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

    def get_skewersx_obj(self, nstride):
        '''
        Return ChollaOnTheFlySkewers_i object of the x-skewers

        Args:
            - n_stride (int): stride cell number between skewers along x-dim
        Return:
            OTFSkewerx (ChollaOnTheFlySkewers_i): skewer object
        '''

        OTFSkewersxHead = ChollaOnTheFlySkewers_iHead(self.nx, self.ny, self.nz,
                                                      nstride, self.xskew_str)

        OTFSkewerx = ChollaOnTheFlySkewers_i(OTFSkewersxHead, self.fPath)

        return OTFSkewerx

