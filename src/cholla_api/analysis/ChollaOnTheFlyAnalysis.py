import numpy as np
import h5py


class ChollaOnTheFlyPowerSpectrumHead:
    '''
    Cholla On The Fly Analysis Power Spectrum Head
    
    Holds information regarding the power spectrum calculation

        Initialized with:
        - dlogk (float): differential step in log k-space
        - nfft (int): number of FFT bins
        - nlos (int): number of line-of-sight cells
        - dvHubble (float): differential Hubble flow velocity

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, dlogk, nfft, nlos, dvHubble):
        self.dlogk = dlogk
        self.n_fft = nfft
        self.n_los = nlos
        self.dvHubble = dvHubble

        self.l_kmin = np.log10( (2 * np.pi) / (self.n_los * self.dvHubble) )
        self.l_kmax = np.log10( (2 * np.pi * (n_fft - 1) ) / (self.n_los * self.dvHubble) )
        self.l_kstart = np.log10(0.99) + self.l_kmin
        self.n_bins = int(1 + ( (self.l_kmax - self.l_kstart) / self.dlogk ) )

    def set_fluxinfo(self, ChollaOTFPowerSpectrum):
        '''
        Set flux attributes for this object

        Args:
            ChollaOTFPowerSpectrum (ChollaOnTheFlyPowerSpectrum): which power spectrum to use to set info
        Returns:
            ...
        '''

        fPath = ChollaOTFPowerSpectrum.fPath
        fObj = h5py.File(fPath, 'r')

        self.HI_flux_mean = fObj['lya_statistics'].attrs[''].item()
        self.HeII_flux_mean = fObj['lya_statistics'].attrs[''].item()
        self.n_skews = fObj['lya_statistics'].attrs[''].item()

        fObj.close()


class ChollaOnTheFlyPhaseSpaceHead:
    '''
    Cholla On The Fly Analysis Phase Space Head

    Holds information regarding the phase space calculation

        Initialized with:
        - density_min (float): minimum density
        - density_max (float): maximum density
        - T_min (float): minimum temperature
        - T_max (float): maximum temperature
        - n_density (int): number of density bins
        - n_temp (int): number of temperature bins

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, density_min, density_max, T_min, T_max, n_density, n_temp):
        self.density_min, self.density_max = density_min, density_max
        self.T_min, self.T_max = T_min, T_max
        self.n_density = n_density
        self.n_temp = n_temp




class ChollaOnTheFlyPowerSpectrum:
    '''
    Cholla On The Fly Analysis Power Spectrum
    
    Holds power spectrum specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFPowerSpectrumHead (ChollaOnTheFlyPowerSpectrumHead): header
            information for power spectrum
        - fPath (str): file path to analysis output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFPowerSpectrumHead, fPath):
       self.OTFPowerSpectrumHead = ChollaOTFPowerSpectrumHead
       self.fPath = fPath
       self.powerspectrum_key = 'p(k)'
       self.kcenters_key = 'k_vals'

    def get_kvals(self, dtype=np.float32):
        '''
        Return the k-centers of the power spectrum

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            kcenters (arr): k mode centers array
        '''

        kcenters = np.zeros(self.OTFPowerSpectrumHead.n_bins, dtype=dtype)
        iter_arr = np.arange(self.OTFPowerSpectrumHead.n_bins, dtype=dtype)
        
        # calculation can be completed in-line, but not as easy to read
        l_kstart = self.OTFPowerSpectrumHead.l_kstart
        dlogk = self.OTFPowerSpectrumHead.dlogk

        kcenters[:] = 10**(l_kstart + (dlogk) * (iter_arr + 1/2) )

        return kcenters 

    def get_powerspectrum(self, dtype=np.float32):
        '''
        Return the power spectrum

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): power spectrum
        '''

        arr = np.zeros(self.OTFPowerSpectrumHead.n_bins, dtype=dtype)
        fObj = h5py.File(self.fPath, 'r')
        arr[:] = fObj['lya_statistics']['power_spectrum'].get(self.powerspectrum_key)
        fObj.close()

        return arr

    def test_kvals(self, dtype=np.float32):
        '''
        Make sure the calculation of k-centers agrees with analysis output

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            (bool): whether a k-value doesn't match up
        '''
        
        kcenters_output = np.zeros(self.OTFPowerSpectrumHead.n_bins, dtype=dtype)
        fObj = h5py.File(self.fPath, 'r')
        arr[:] = fObj['lya_statistics']['power_spectrum'].get(self.kcenters_key)
        fObj.close()

        # TODO: MAKE SURE THE TWO ARRAYS ARE OF TEH SAME SHAPE TOO!
        kcenters = self.get_kvals(dtype=dtype)
        if (not np.array_equal(kcenters_output, kcenters) ):
            # find the indices where they don't match up
            failingk_inds = np.argwhere(kcenters_output != kcenters).flatten()

            print("Failed k-value calculation")
            for k_ind in failingk_inds:
                print(f"\tFailed at index: {k_ind:.0f}")
                print(f"\t\tCalculated k[{k_ind:.0f}] = {kcenters[k_ind]:.5e}")
                print(f"\t\tExpected k[{k_ind:.0f}] = {kcenters_output[k_ind]:.5e}")
            return False

        return True


class ChollaOnTheFlyPhaseSpace:
    '''
    Cholla On The Fly Analysis Phase Space
    
    Holds phase space specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFPhaseSpaceHead (ChollaOnTheFlyPhaseSpaceHead): header
            information associated with phase space
        - fPath (str): file path to analysis output

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, ChollaOTFPhaseSpaceHead, fPath):
       self.OTFPhaseSpaceHead = ChollaOTFPhaseSpaceHead
       self.fPath = fPath

    def get_phasespace(self, dtype=np.float32):
        '''
        Return the phase space array

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): phase space
        '''

        # can be computed in-line, but harder to read
        phase_shape = (self.OTFPhaseSpaceHead.n_density, 
                       self.OTFPhaseSpaceHead.n_temp)
        arr = np.zeros(phase_shape, dtype=dtype)
        fObj = h5py.File(self.fPath, 'r')
        arr[:] = fObj['phase_diagram'].get('data')
        fObj.close()

        return arr
    

class ChollaOnTheFlyAnalysis:
    '''
    Cholla On The Fly Analysis
    
    Holds on-the-fly analysis specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFAnalysisHead (ChollaOnTheAnalysisHead): header information 
            for analysis output
        - AnalysisPath (str): directory path to analysis output files
        - ChollaGrid (ChollaGrid): grid holding domain information

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, nAnalysis, AnalysisPath, ChollaGrid):
        self.OTFAnalysisfPath = AnalysisPath + '/' + str(nAnalysis) + '_analysis.h5'

        # current implementation only works for cube grid
        assert ChollaGrid.nx_global == ChollaGrid.ny_global
        assert ChollaGrid.nx_global == ChollaGrid.nz_global
        assert ChollaGrid.Lx == ChollaGrid.Ly
        assert ChollaGrid.Lx == ChollaGrid.Lz

        self.n_fft = int(ChollaGrid.nx_global / 2 + 1)
        self.n_los = ChollaGrid.nx_global
        dx = ChollaGrid.dx / 1e3 # convert kpc --> Mpc

        # set cosmology params
        self.set_cosmoinfo()

        # grab current hubble param & info needed to calculate hubble flow
        H = self.get_currH()
        cosmoh = self.H0 / 100
        dxproper = dx * self.current_a / cosmoh

        # calculate hubble flow through a cell
        self.dvHubble = H * dxproper

        # hard coded in Cholla: src/analysis/phase_diagram.cpp
        self.n_dens   = 1000
        self.n_temp   = 1000
        self.dens_min = 1e-3
        self.dens_max = 1e6
        self.temp_min = 1e0
        self.temp_max = 1e8


    def set_cosmoinfo(self):
        '''
        Set cosmological attributes for this object

        Args:
            ChollaOTFAnalysis (ChollaOnTheFlyAnalysis): specific Analysis object to use
        Returns:
            ...
        '''

        fObj = h5py.File(self.OTFAnalysisfPath, 'r')

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
        DE_factor *= np.exp(-3. * (self.wa * (1. - self.current_a)))
        
        H0_factor = (self.Omega_R / a4) + (self.Omega_M / a3) 
        H0_factor += (self.Omega_K / a2) + (self.Omega_L * DE_factor)
       
        return self.H0 * H0_factor

    def get_powerspectrum_obj(self, dlogk):
        '''
        Return ChollaOnTheFlyPowerSpectrum object of this analysis

        Args:
            dlogk (float): differential step in log k-space
        Return:
            OTFPowerSpectrum (ChollaOnTheFlyPowerSpectrum): power spectrum object
        '''

        OTFPowerSpectrumHead = ChollaOnTheFlyPowerSpectrumHead(dlogk, self.n_fft, 
                                                               self.n_los, self.dvHubble)
        
        OTFPowerSpectrum = ChollaOnTheFlyPowerSpectrum(OTFPowerSpectrumHead, 
                                                       self.OTFAnalysisfPath)
        # place flux info
        OTFPowerSpectrumHead.set_fluxinfo(OTFPowerSpectrum)

        return OTFPowerSpectrum

    def get_phasespace_obj(self):
        '''
        Return ChollaOnTheFlyPhaseSpace object of this analysis

        Args:
            ...
        Return:
            OTFPhaseSpace (ChollaOnTheFlyPhaseSpace): phase space object            
        '''

        OTFPhaseSpaceHead = ChollaOnTheFlyPhaseSpaceHead(self.dens_min, self.dens_max, 
                                                         self.temp_min, self.temp_max, 
                                                         self.n_dens, self.n_temp)

        OTFPhaseSpace = ChollaOnTheFlyPhaseSpace(OTFPhaseSpaceHead, self.OTFAnalysisfPath)

        return OTFPhaseSpace

    def get_phasespace(self, dtype=np.float32):
        '''
        Return the 2D phase space array

        Args:
            dtype (np type): (optional) numpy precision to use
        Return:
            (arr): phase space array
        '''

        return self.get_phasespace_obj().get_phasespace(dtype=dtype)
        
    def get_kvals(self, dlogk, dtype=np.float32):
        '''
        Return the k values used in the power spectrum

        Args:
            dlogk (float): differential step in log k-space
            dtype (np type): (optional) numpy precision to use
        Return:
            (arr): k-mode values
        '''

        return self.get_powerspectrum_obj(dlogk).get_kvals(dtype=dtype)

    def get_powerspectrum(self, dlogk, dtype=np.float32):
        '''
        Return the power spectrum

        Args:
            dlogk (float): differential step in log k-space
            dtype (np type): (optional) numpy precision to use
        Return:
            (arr): power spectrum
        '''

        return self.get_powerspectrum_obj(dlogk).get_powerspectrum(dtype=dtype)





