import numpy as np
from scipy.special import erf

from cholla_api.analysis.ChollaHydroCalculator import ChollaHydroCalculator
from cholla_api.analysis.ChollaCosmoCalculator import ChollaSnapCosmologyHead, ChollaCosmoCalculator


class ChollaSkewerCosmoCalculator:
    '''
    Cholla Skewer Calculator object
        Serves as a specific implementaiton of a Cholla Cosmological Calculator
            for a skewer.

        Initialized with:
            scale_factor (float): scale factor
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units
            n_los (int): number of cells along line-of-sight
            dx (float): comoving distance between cells (kpc)
            dtype (np type): (optional) numpy precision to initialize output arrays
        
        Objects including ghost cells are suffixed with _ghost

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, scale_factor, cosmoHead, n_los, dx, dtype=np.float32):
        self.n_los = n_los
        self.n_ghost = int(0.1 * n_los) # take 10% from bruno
        self.dx = dx
        self.a = scale_factor

        # number of line-of-sight cells including ghost cells
        self.n_los_ghost = self.n_los + 2 * self.n_ghost

        # create ChollaCosmoCalc object
        self.snapCosmoHead = ChollaSnapCosmologyHead(self.a, cosmoHead)
        calc_dims, calc_dims_ghost = (self.n_los,), (self.n_los_ghost,)
        self.snapCosmoCalc = ChollaCosmoCalculator(self.snapCosmoHead, calc_dims, dtype=dtype)
        self.snapCosmoCalc_ghost = ChollaCosmoCalculator(self.snapCosmoHead, calc_dims_ghost, dtype=dtype)

        # create HydroCalc objects
        self.hydroCalc = ChollaHydroCalculator(calc_dims, dtype=dtype)
        self.hydroCalc_ghost = ChollaHydroCalculator(calc_dims_ghost, dtype=dtype)

        # calculate Hubble flow through one cell
        dvHubble = self.snapCosmoHead.dvHubble(self.dx) # [km s-1]
        self.dvHubble_cgs = dvHubble * self.snapCosmoHead.cosmoHead.km_cgs # [cm s-1]

        # create Hubble flow arrays along left, right, and center of each cell
        # prepend and append ghost cells
        self.vHubbleL_ghost_cgs = np.arange(-self.n_ghost, self.n_ghost + self.n_los) * self.dvHubble_cgs
        self.vHubbleR_ghost_cgs = self.vHubbleL_ghost_cgs + self.dvHubble_cgs
        self.vHubbleC_ghost_cgs = self.vHubbleL_ghost_cgs + 0.5 * self.dvHubble_cgs


    def extend_ghostcells(self, arr):
        '''
        Extend an array with ghost cels, enforcing periodic boundary conditions
        
        Args:
            arr (arr): array to extend
        Returns:
            arr_ghost (arr): extended array
        '''

        # create array with correct size, inherit arr datatype
        arr_ghost = np.zeros(self.n_los_ghost, dtype=arr.dtype)
    
        # copy over the real data
        arr_ghost[self.n_ghost : self.n_ghost + self.n_los] = arr[:]
        
        # extend far edge of arr to near edge of ghost array
        arr_ghost[ : self.n_ghost] = arr[-self.n_ghost : ]

        # extend near edge of arr to far edge of ghost array
        arr_ghost[-self.n_ghost : ] = arr[ : self.n_ghost]

        return arr_ghost


    def optical_depth_Hydrogen(self, densityHI, velocity_pec, temp, num_sigs=0):
        '''
        Compute the optical depth for each cell along the line-of-sight

        Args:
            densityHI (arr): ionized Hydrogen comoving density [h2 Msun kpc-3]
            velocity_pec (arr): peculiar velocity [km s-1]
            temp (arr): temperature [K]
            num_sigs (int): (optional) number of standard deviations to around mean
        Returns:
            tau (arr): optical depth for each cell
        '''
        assert densityHI.size == self.n_los
        assert velocity_pec.size == self.n_los
        assert temp.size == self.n_los

        # convert comoving density to physical density then to cgs
        densityHI_phys = self.snapCosmoCalc.physical_density(densityHI)
        densityHI_phys_cgs = self.snapCosmoCalc.density_cosmo2cgs(densityHI_phys) # [g cm-3]

        # calculate column number density & extend to ghost cells
        nHI_phys_cgs = densityHI_phys_cgs / self.hydroCalc.mp # [cm-3]
        nHI_phys_ghost_cgs = self.extend_ghostcells(nHI_phys_cgs)

        # convert peculiar velocity to cgs values & extend to ghost cells
        velocity_pec_cgs = self.snapCosmoCalc.velocity_cosmo2cgs(velocity_pec)
        velocity_pec_ghost_cgs = self.extend_ghostcells(velocity_pec_cgs)
        # convert peculiar to physical velocity by adding Hubble flow
        velocity_phys_ghost_cgs = velocity_pec_ghost_cgs + self.vHubbleC_ghost_cgs # [cm s-1]

        # calculate doppler broadening param & extend to ghost cells
        doppler_param_cgs = self.hydroCalc.Doppler_param_Hydrogen(temp) # [cm s-1]
        doppler_param_ghost_cgs = self.extend_ghostcells(doppler_param_cgs)

        # calculate Ly-alpha interaction cross section
        sigma_Lya = np.pi * self.hydroCalc.e * self.hydroCalc.e / self.hydroCalc.me # [cm3 g1 s-2 / g] = [cm3 s-2]
        sigma_Lya = sigma_Lya * self.hydroCalc.lambda_Lya / self.hydroCalc.c # [cm3 s-2 * cm / (cm s-1)] = [cm3 s-1]
        sigma_Lya = sigma_Lya / self.snapCosmoHead.Hubble_cgs # [cm3 s-1 / (s-1)] = [cm3]
        f_12 = 0.416 # oscillator strength
        sigma_Lya *= f_12

        # initialize optical depths
        tau_ghost = self.snapCosmoCalc_ghost.create_arr()
        

        if num_sigs == 0:
            # OLD IMPLEMENTATION
            for losid in range(self.n_los_ghost):
                vH_L, vH_R = self.vHubbleL_ghost_cgs[losid], self.vHubbleR_ghost_cgs[losid]
                # calculate line center shift in terms of broadening scale
                y_L = (vH_L - velocity_phys_ghost_cgs) / doppler_param_ghost_cgs
                y_R = (vH_R - velocity_phys_ghost_cgs) / doppler_param_ghost_cgs
                # [cm3 * # density] = [cm3 * cm-3] = []
                tau_ghost[losid] = sigma_Lya * np.sum(nHI_phys_ghost_cgs * (erf(y_R) - erf(y_L))) / 2.0
        else:        
            # NEW IMPLEMENTATION -- dynamic scaling of 5 sigs
            five_sigs = (num_sigs) * (doppler_param_ghost_cgs)
            
            vHC_fivesig_upp_all = self.vHubbleC_ghost_cgs + five_sigs
            vHC_fivesig_low_all = self.vHubbleC_ghost_cgs - five_sigs

            for losid in range(self.n_los_ghost):
                vH_L, vH_R = self.vHubbleL_ghost_cgs[losid], self.vHubbleR_ghost_cgs[losid]
                # evaluate the five-sig mask here
                vHC_fivesig_upp = vHC_fivesig_upp_all[losid]
                vHC_fivesig_low = vHC_fivesig_low_all[losid]
                fivesig_mask = (velocity_phys_ghost_cgs < vHC_fivesig_upp) & (velocity_phys_ghost_cgs > vHC_fivesig_low)
                vel_fivesig = velocity_phys_ghost_cgs[fivesig_mask]
                doppler_fivesig = doppler_param_ghost_cgs[fivesig_mask]
                nHI_fivesig = nHI_phys_ghost_cgs[fivesig_mask]
                # calculate line center shift in terms of broadening scale
                y_L = (vH_L - vel_fivesig) / doppler_fivesig
                y_R = (vH_R - vel_fivesig) / doppler_fivesig
                tau_ghost[losid] = (sigma_Lya / 2.0) * np.sum(nHI_fivesig * (erf(y_R) - erf(y_L)))


        # clip edges
        tau = tau_ghost[self.n_ghost : -self.n_ghost]

        return tau


