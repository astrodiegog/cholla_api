import numpy as np
import h5py

from cholla_api.analysis.ChollaFluxPowerSpectrum import ChollaFluxPowerSpectrumHead

class ChollaOnTheFlyFluxPowerSpectrum:
    '''
    Cholla Flux Power Spectrum
    
    Holds information to perform the power spectrum calculation from the skewers
        saved on-the-fly.

        Initialized with:
        - dlogk (float): differential step in log k-space
        - chOTFSkewers (ChollaOnTheFlySkewers): interacts with skewer output files 

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, dlogk, chOTFSkewers):
        self.OTFSkewers = chOTFSkewers
        self.dlogk = dlogk

    def get_alllocalopticaldepth(self, OTFSkewers_i, dtype=np.float32):
        '''
        Return the local optical depth dataset for all skewers.

        Args:
            OTFSkewers_i (ChollaOnTheFlySkewers_i): skewer object
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        fObj = h5py.File(OTFSkewers_i.fPath, 'r')

        # ensure that local optical depth is a dataset 
        local_opticaldepth_key = 'taucalc_local'
        assert local_opticaldepth_key in fObj[OTFSkewers_i.OTFSkewersiHead.skew_key].keys()

        arr = np.zeros((OTFSkewers_i.OTFSkewersiHead.n_skews, OTFSkewers_i.OTFSkewersiHead.n_i), dtype=dtype)
        arr[:,:] = fObj[OTFSkewers_i.OTFSkewersiHead.skew_key].get(local_opticaldepth_key)[:, :]
        fObj.close()

        return arr

    def get_FPS_x(self, precision=np.float64):
        '''
        Return the Flux Power Spectrum along the x-axis

        Args:
            precision (np type): (optional) numpy precision to use
        Return:
            (arr): k mode edges array
            (arr): mean transmitted flux power spectrum within kmode edges
        '''
        # grab x-skewer object
        OTFSkewers_x = self.OTFSkewers.get_skewersx_obj()

        # create Flux Power Spectrum Head object
        FPSHead_x = ChollaFluxPowerSpectrumHead(self.dlogk, OTFSkewers_x.OTFSkewersiHead.n_i,
                                                self.OTFSkewers.dvHubble_x)

        # grab local optical depths
        local_opticaldepth = self.get_alllocalopticaldepth(OTFSkewers_x, precision)

        # return flux power spectrum along x-axis
        return FPSHead_x.get_FPS(local_opticaldepth, precision)

    def get_FPS_y(self, precision=np.float64):
        '''
        Return the Flux Power Spectrum along the y-axis

        Args:
            precision (np type): (optional) numpy precision to use
        Return:
            (arr): k mode edges array
            (arr): mean transmitted flux power spectrum within kmode edges
        '''
        # grab y-skewer object
        OTFSkewers_y = self.OTFSkewers.get_skewersy_obj()

        # create Flux Power Spectrum Head object
        FPSHead_y = ChollaFluxPowerSpectrumHead(self.dlogk, OTFSkewers_y.OTFSkewersiHead.n_i,
                                                self.OTFSkewers.dvHubble_y)

        # grab local optical depths
        local_opticaldepth = self.get_alllocalopticaldepth(OTFSkewers_y, precision)

        # return flux power spectrum along y-axis
        return FPSHead_y.get_FPS(local_opticaldepth, precision)


    def get_FPS_z(self, precision=np.float64):
        '''
        Return the Flux Power Spectrum along the z-axis

        Args:
            precision (np type): (optional) numpy precision to use
        Return:
            (arr): k mode edges array
            (arr): mean transmitted flux power spectrum within kmode edges
        '''
        # grab z-skewer object
        OTFSkewers_z = self.OTFSkewers.get_skewersz_obj()

        # create Flux Power Spectrum Head object
        FPSHead_z = ChollaFluxPowerSpectrumHead(self.dlogk, OTFSkewers_z.OTFSkewersiHead.n_i,
                                                self.OTFSkewers.dvHubble_z)

        # grab local optical depths
        local_opticaldepth = self.get_alllocalopticaldepth(OTFSkewers_z, precision)

        # return flux power spectrum along z-axis
        return FPSHead_z.get_FPS(local_opticaldepth, precision)
