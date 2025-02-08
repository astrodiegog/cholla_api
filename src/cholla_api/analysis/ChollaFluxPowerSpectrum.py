import numpy as np

class ChollaFluxPowerSpectrumHead:
    '''
    Cholla Flux Power Spectrum Head
    
    Holds information regarding the power spectrum calculation

        Initialized with:
        - nlos (int): number of line-of-sight cells
        - dvHubble (float): differential Hubble flow velocity across a cell

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, nlos, dvHubble):
        self.n_los = nlos
        self.n_fft = int(self.n_los / 2 + 1)
        self.dvHubble = dvHubble

        # calculate Hubble flow across entire box (max velocity)
        self.u_max = self.dvHubble * self.n_los

        self.l_kmin = np.log10( (2. * np.pi) / (self.u_max) )
        self.l_kmax = np.log10( (2. * np.pi * (self.n_fft - 1.) ) / (self.u_max) )

    def get_kvals_fft(self, dtype=np.float32):
        '''
        Return k-modes from the Fourier Transform

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            kcenters_fft (arr): k mode centers array
        '''

        kcenters_fft = np.zeros(self.n_fft, dtype=dtype)
        iter_arr = np.arange(self.n_fft, dtype=dtype)

        kcenters_fft[:] = (2. * np.pi * iter_arr) / (self.u_max)

        return kcenters_fft


    def get_FPS(self, local_opticaldepths, precision=np.float64):
        '''
        Return the Flux Power Spectrum given the local optical depths.
            Expect 2-D array of shape (number skewers, line-of-sight cells)

        Args:
            local_opticaldepths (arr): local optical depths of all skewers
            mean_flux (float): mean flux to scale deviations
            precision (np type): (optional) numpy precision to use
        Return:
            kmode_fft (arr): Fourier Transform k mode array
            P_k_mean (arr): mean transmitted flux power spectrum within kmode edges
        '''
        assert local_opticaldepths.ndim == 2
        assert local_opticaldepths.shape[1] == self.n_los

        n_skews = local_opticaldepths.shape[0]

        # calculate local transmitted flux (& its mean)
        fluxes = np.exp(-local_opticaldepths)
        flux_mean = np.mean(fluxes)

        # initialize total power array & temporary FFT array
        delta_F_avg = np.zeros(self.n_fft, dtype=precision)
        P_k_tot = np.zeros(self.n_fft, dtype=precision)

        for nSkewerID in range(n_skews):
            # clean out temporary FFT array
            delta_F_avg[:] = 0.

            # calculate flux fluctuation 
            dFlux_skew = (fluxes[nSkewerID] - flux_mean) / flux_mean

            # perform fft & calculate amplitude of fft
            fft = np.fft.rfft(dFlux_skew)
            fft2 = (fft.imag * fft.imag) + (fft.real * fft.real)

            # take avg & scale by umax
            delta_F_avg = fft2 / self.n_los / self.n_los
            P_k = self.u_max * delta_F_avg
            P_k_tot += P_k

        # average out by the number of skewers
        P_k_mean = P_k_tot / n_skews

        # grab k-mode bin edges
        kmode_fft = self.get_kvals_fft(precision)

        return (kmode_fft, P_k_mean)


    def rebin_FPS(self, dlogk, dtype_bin=np.uint64, precision=np.float32):
        '''
        Bin the Flux Power Spectrum along some differential log-space scaling

        Args:
            dlogk (float): differential step in log k-space
            dtype_bin (np type): (optional) numpy precision to use for bin ids
            precision (np type): (optional) numpy precision to use for k values
        Returns:
            (arr): k-mode edges array
            (arr): number of FFT bins falling in each k-mode bin
        '''

        l_kstart = np.log10(0.99) + self.l_kmin
        n_bins = int(1 + ( (self.l_kmax - l_kstart) / dlogk ) )

        # grab fft kvalues
        kvals_fft = self.get_kvals_fft(dtype=precision)

        # find the indices that describe where the k-mode FFT bins fall within kval_edges (from dlogk)
        fft_binids = np.zeros(self.n_fft, dtype=dtype_bin)
        fft_binids_float = (np.log10(kvals_fft) - l_kstart) / dlogk
        fft_binids[:] = np.floor(fft_binids_float)

        # count the number of fft modes that fall within each dlogk bin id
        kedges_nfft_bins = np.zeros(n_bins, dtype=dtype_bin)
        for fft_binid in fft_binids[1:]:
            kedges_nfft_bins[fft_binid] += 1
        
        # calculate new kedges
        kedges = np.zeros(n_bins + 1, dtype=precision)
        iter_arr = np.arange(n_bins + 1, dtype=precision)
        kedges[:] = 10**(l_kstart + (dlogk * iter_arr) )

        return kedges, fft_binids

