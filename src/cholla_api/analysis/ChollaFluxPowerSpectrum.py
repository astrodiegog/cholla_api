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


