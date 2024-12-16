import numpy as np

class ChollaFluxPowerSpectrumHead:
    '''
    Cholla Flux Power Spectrum Head
    
    Holds information regarding the power spectrum calculation

        Initialized with:
        - dlogk (float): differential step in log k-space
        - nlos (int): number of line-of-sight cells
        - dvHubble (float): differential Hubble flow velocity across a cell

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, dlogk, nlos, dvHubble):
        self.dlogk = dlogk
        self.n_los = nlos
        self.n_fft = int(self.n_los / 2 + 1)
        self.dvHubble = dvHubble

        # calculate Hubble flow across entire box (max velocity)
        self.u_max = self.dvHubble * self.n_los

        self.l_kmin = np.log10( (2. * np.pi) / (self.u_max) )
        self.l_kmax = np.log10( (2. * np.pi * (self.n_fft - 1.) ) / (self.u_max) )
        self.l_kstart = np.log10(0.99) + self.l_kmin
        self.n_bins = int(1 + ( (self.l_kmax - self.l_kstart) / self.dlogk ) )


    def get_kvals(self, dtype=np.float32):
        '''
        Return the k-centers of the power spectrum

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            kcenters (arr): k mode centers array
        '''

        kcenters = np.zeros(self.n_bins, dtype=dtype)
        iter_arr = np.arange(self.n_bins, dtype=dtype)

        kcenters[:] = 10**(self.l_kstart + (self.dlogk) * (iter_arr + 0.5) )

        return kcenters

    def get_kvals_edges(self, dtype=np.float32):
        '''
        Return the k-edges of the power spectrum

        Args:
            dtype (np type): (optional) numpy precision to use
        Returns:
            kedges (arr): k mode edges array
        '''

        kedges = np.zeros(self.n_bins + 1, dtype=dtype)
        iter_arr = np.arange(self.n_bins + 1, dtype=dtype)

        kedges[:] = 10**(self.l_kstart + (self.dlogk * iter_arr) )

        return kedges

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

    def get_fft_binids(self, dtype_bin=np.int64, dtype_calc=np.float32, useforloop=True):
        '''
        Return the indices that the k-mode fft bins land on within kvals_edges

        Args:
            dtype_bin (np type): (optional) numpy precision to use for returned array
            dtype_calc (np type): (optional) numpy precision to use for calculations
            useforloop (bool): (optional) whether to use for-loop or not  
        Returns:
            fft_binids (arr): indices where fft k-mode lands wrt kvals edges
        '''

        fft_binids = np.zeros(self.n_fft, dtype=dtype_bin)

        # grab fft kvalues
        kvals_fft = self.get_kvals_fft(dtype=dtype_calc)

        if useforloop:
            # grab edges for comparison
            kvals_edges = self.get_kvals_edges(dtype=dtype_calc)
            for bin_id_fft in range(self.n_fft):
                edge_greater_fft = np.argwhere(kvals_fft[bin_id_fft] < kvals_edges).flatten()
                # okay to flatten bc we know kvals_fft and kvals_edges are 1D arrays
                if edge_greater_fft.size > 0 :
                    # ensure we're indexing into a non-empty array
                    fft_binids[bin_id_fft] = edge_greater_fft[0] - 1
        else:
            fft_binids_float = (np.log10(kvals_fft) - self.l_kstart) / self.dlogk
            fft_binids[:] = np.floor(fft_binids_float)

        return fft_binids


    def get_FPS(self, local_opticaldepths, precision=np.float64):
        '''
        Return the Flux Power Spectrum given the local optical depths.
            Expect 2-D array of shape (number skewers, line-of-sight cells)

        Args:
            local_opticaldepths (arr): local optical depths of all skewers
            precision (np type): (optional) numpy precision to use
        Return:
            kmode_edges (arr): k mode edges array
            P_k_mean (arr): mean transmitted flux power spectrum within kmode edges
        '''
        assert local_opticaldepths.ndim == 2
        assert local_opticaldepths.shape[1] == self.n_los

        n_skews = local_opticaldepths.shape[0]

        # find the indices that describe where the k-mode FFT bins fall within kval_edges (from dlogk)
        fft_binids = self.get_fft_binids(dtype_bin=np.int64, dtype_calc=np.float64)

        # find number of fft modes that fall in requested dlogk bin id (used later to average total power in dlogk bin)
        hist_n = np.zeros(self.n_bins, dtype=precision)
        for bin_id in fft_binids[1:]:
            hist_n[bin_id] += 1.
        # (protect against dividing by zero)
        hist_n[hist_n == 0] = 1.

        # calculate local transmitted flux (& its mean)
        fluxes = np.exp(-local_opticaldepths)
        flux_mean = np.mean(fluxes)

        # initialize total power array & temporary FFT array
        hist_PS_vals = np.zeros(self.n_bins, dtype=precision)
        P_k_tot = np.zeros(self.n_bins, dtype=precision)

        for nSkewerID in range(n_skews):
            # clean out temporary FFT array
            hist_PS_vals[:] = 0.

            # calculate flux fluctuation 
            dFlux_skew = fluxes[nSkewerID] / flux_mean

            # perform fft & calculate amplitude of fft
            fft = np.fft.rfft(dFlux_skew)
            fft2 = (fft.imag * fft.imag + fft.real * fft.real) / self.n_los / self.n_los

            # add power for each fft mode
            hist_PS_vals[fft_binids[1:]] += fft2[1:]

            # take avg & scale by umax
            delta_F_avg = hist_PS_vals / hist_n
            P_k = self.u_max * delta_F_avg
            P_k_tot += P_k

        # average out by the number of skewers
        P_k_mean = P_k_tot / n_skews

        # grab k-mode bin edges
        kmode_edges = self.get_kvals_edges(precision)

        return (kmode_edges, P_k_mean)





