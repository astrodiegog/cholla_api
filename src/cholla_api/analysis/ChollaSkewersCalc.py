import numpy as np
from scipy.special import erf

from cholla_api.analysis.ChollaCosmoGridCalculator import ChollaCosmoCalculator


class ChollaSkewerCosmoCalculator:
    '''
    Cholla Skewer Calculator object
        Serves as a specific implementaiton of a Cholla Cosmological Calculator
            for a skewer.

        Initialized with:
            snapHead (ChollaSnapHead): provides current redshift
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units
            n_los (int): number of cells along line-of-sight
            dx (float): comoving distance between cells (kpc)
            dtype (np type): (optional) numpy precision to initialize output arrays
        
        Objects including ghost cells are suffixed with _ghost

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, snapHead, cosmoHead, n_los, dx, dtype=np.float32):
        self.n_los = n_los
        self.n_ghost = int(0.1 * n_los) # take 10% from bruno
        self.dx = dx
        self.dtype = dtype

        self.snapHead = snapHead

        # number of line-of-sight cells including ghost cells
        self.n_los_ghost = self.n_los + 2 * self.n_ghost

        self.cosmoCalc_ghost = ChollaCosmoCalculator(self.snapHead, cosmoHead)

        # calculate Hubble flow through one cell
        dvHubble = self.cosmoCalc_ghost.dvHubble(self.dx) # [km s-1]
        self.dvHubble_cgs = dvHubble * self.cosmoCalc_ghost.cosmoHead.km_cgs # [cm s-1]

        # create Hubble flow arrays along left, right, and center of each cell
        # prepend and append ghost cells
        #self.vHubbleL_ghost = np.arange(-self.n_ghost, self.n_ghost + self.n_los) * self.dvHubble
        #self.vHubbleR_ghost = vHubbleL + self.dvHubble

        # create Hubble flow array at center of each cell
        self.vHubbleCenter_ghost_cgs = (np.arange(-self.n_ghost, self.n_ghost + self.n_los) + 0.5) * self.dvHubble_cgs # [cm s-1]

        # calculate & attach current Hubble rate in [km s-1 Mpc-1] and [s-1]
        self.Hubble = self.cosmoCalc_ghost.Hubble()
        self.Hubble_cgs = self.Hubble * self.cosmoCalc_ghost.cosmoHead.km_cgs / self.cosmoCalc_ghost.cosmoHead.Mpc_cgs # in cgs [s-1]


    def extend_ghostcells(self, arr):
        '''
        Extend an array with ghost cels, enforcing periodic boundary conditions
        
        Args:
            arr (arr): array to extend
        Returns:
            arr_ghost (arr): extended array
        '''

        # create array with correct size
        arr_ghost = self.cosmoCalc_ghost.create_arr()
    
        # copy over the real data
        arr_ghost[self.n_ghost : self.n_ghost + self.n_los] = arr[:]
        
        # extend far edge of arr to near edge of ghost array
        arr_ghost[ : self.n_ghost] = arr[-self.n_ghost : ]

        # extend near edge of arr to far edge of ghost array
        arr_ghost[-self.n_ghost : ] = arr[ : self.n_ghost]

        return arr_ghost


    def optical_depth_Hydrogen(self, densityHI, velocity_pec, temp):
        '''
        Compute the optical depth for each cell along the line-of-sight

        Args:
            densityHI (arr): ionized Hydrogen comoving density [h2 Msun kpc-3]
            velocity_pec (arr): peculiar velocity [km s-1]
            temp (arr): temperature [K]
        Returns:
            tau (arr): optical depth for each cell
        '''
        assert densityHI.size == self.n_los
        assert velocity_pec.size == self.n_los
        assert temp.size == self.n_los

        # introduce hydro calculator
        hydroCalc = ChollaHydroCalculator(self.n_los, dtype=self.cosmoCalc_ghost.dtype)
        
        # convert comoving density to phyiscal density
        densityHI_phys = self.cosmoCalc_ghost.physical_density(densityHI) # [h2 Msun kpc-3]
        # convert density to cgs
        h_cosmo2 = self.cosmoCalc_ghost.cosmoHead.h_cosmo * self.cosmoCalc_ghost.cosmoHead.h_cosmo
        densityHI_phys_cgs = densityHI_phys * self.cosmoCalc_ghost.cosmoHead.Msun_cgs * h_cosmo2 / self.cosmoCalc_ghost.cosmoHead.kpc3_cgs # [g cm-3]
        # calculate column number density
        nHI_phys_cgs = densityHI_phys_cgs / hydroCalc.mp # [cm-3]
        # extend number density
        nHI_phys_ghost_cgs = self.extend_ghostcells(nHI_phys_cgs)


        # extend peculiar velocity
        velocity_pec_ghost = self.extend_ghostcells(velocity_pec)
        # convert velocity to cgs
        velocity_pec_ghost_cgs = velocity_pec_ghost * self.cosmoCalc_ghost.cosmoHead.km_cgs # [cm s-1]
        # convert to physical velocity by adding Hubble flow
        velocity_phys_ghost_cgs = velocity_pec_ghost + self.vHubbleCenter_ghost_cgs # [cm s-1]
        

        # calculate doppler broadening param
        hydroCalc = ChollaHydroCalculator(self.n_los, dtype=self.cosmoCalc_ghost.dtype)
        doppler_param_cgs = hydroCalc.Doppler_param_Hydrogen(temp) # [cm s-1]
        # extend doppler param
        doppler_param_ghost_cgs = self.extend_ghostcells(doppler_param_cgs)


        # calculate Ly-alpha interaction cross section
        sigma_Lya = np.pi * hydroCalc.e * hydroCalc.e / hydroCalc.me # [cm3 g1 s-2 / g] = [cm3 s-2]
        sigma_Lya = sigma_Lya * hydroCalc.lambda_Lya / hydroCalc.c # [cm3 s-2 * cm / (cm s-1)] = [cm3 s-1]
        sigma_Lya = sigma_Lya / self.Hubble_cgs # [cm3 s-1 / (s-1)] = [cm3]
        f_12 = 0.416 # oscillator strength
        sigma_Lya *= f_12

        # initialize optical depths
        tau_ghost = np.zeros(self.n_los_ghost, dtype=self.dtype)
        
        # fill optical depths with Ly-alpha x-section
        tau_ghost += sigma_Lya # [cm3]

        for losid in range(self.self.n_los_ghost):
            # Hubble flow velocity at center of cell
            vH_losid = self.vHubbleCenter_ghost_cgs[los_id]
            # find Hubble flow velocity at left and right interfaces of cell
            vH_L, vH_R = vH_losid - 0.5 * self.dvHubble_cgs, vH_losid + 0.5 * self.dvHubble_cgs
            # calculate line center shift in terms of broadening scale
            y_L = (vH_L - velocity_phys_ghost_cgs) / doppler_param_ghost_cgs
            y_R = (vH_R - velocity_phys_ghost_cgs) / doppler_param_ghost_cgs
            # [cm3 * # density] = [cm3 * cm-3] = []
            tau_ghost[losid] *= np.sum(nHI_phys_ghost_cgs * erf(y_R) - erf(y_R)) / 2.0

        # clip edges
        tau = tau_ghost[n_ghost:-n_ghost]

        return tau


