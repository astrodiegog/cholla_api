import numpy as np

class ChollaCosmologyHead:
    '''
    Cholla Cosmology Head
        Serves as a header object that holds information that helps define a
            specific cosmology
        
        Initialized with:
        - OmegaM (float): present-day energy density parameter for matter
        - OmegaR (float): present-day energy density parameter for radiation
        - OmegaK (float): present-day energy density parameter for spatial curvature
        - OmegaL (float): present-day energy density parameter for dark energy
        - w0 (float): constant term in dark energy equation of state
        - wa (float): linear term in dark energy equation of state
        - H0 (float): present-day Hubble parameter in units of [km / s / Mpc]

    '''
    
    def __init__(self, OmegaM, OmegaR, OmegaK, OmegaL, w0, wa, H0):
        
        # start with constants !
        self.Msun_cgs = 1.98847e33 # Solar Mass in grams
        self.kpc_cgs = 3.0857e21 # kiloparsecs in centimeters
        self.Mpc_cgs = self.kpc_cgs * 1.e3 # Megaparsecs in centimeters
        self.km_cgs = 1.e5 # kilometers in centimeters
        self.kyr_cgs = 3.15569e10 # kilo-years in seconds
        self.Myr_cgs = self.kyr_cgs * 1.e3 # mega-years in seconds
        self.Gyr_cgs = self.Myr_cgs * 1.e3 # giga-years in seconds

        self.G_cgs = 6.67259e-8 # gravitational constant in cgs [cm3 g-1 s-2]
        self.G_cosmo = self.G_cgs / self.km_cgs / self.km_cgs / self.kpc_cgs * self.Msun_cgs # gravitational constant in cosmological units [kpc (km2 s-2) Msun-1]

        self.kpc3_cgs = self.kpc_cgs * self.kpc_cgs * self.kpc_cgs
        self.Mpc3_cgs = self.Mpc_cgs * self.Mpc_cgs * self.Mpc_cgs

        # present-day energy density for matter, radiation, curvature, and Dark Energy
        self.OmegaM = OmegaM
        self.OmegaR = OmegaR
        self.OmegaK = OmegaK
        self.OmegaL = OmegaL

        # Dark Energy equation of state like w(a) = w0 + wa(1-a)
        self.w0, self.wa = w0, wa

        # present-day hubble parameter
        self.H0 = H0 # in [km s-1 Mpc-1]
        self.H0_cgs = self.H0 * self.km_cgs / self.Mpc_cgs # in cgs [s-1]
        self.H0_cosmo = self.H0 / 1.e3 # in cosmological units [km s-1 kpc-1]

        # dimensionless hubble parameter
        self.h_cosmo = self.H0 / 100.

        # Hubble time (1/H0)
        self.t_H0_cgs = 1. / self.H0_cgs # in seconds
        self.t_H0_gyrs = self.t_H0_cgs / self.Gyr_cgs # in Gyrs
        self.t_H0_cosmo  = self.t_H0_cgs * self.km_cgs / self.kpc_cgs # in cosmological units [s kpc km-1]

        # critical density in units of [g cm-3]
        self.rho_crit0_cgs = 3. * self.H0_cgs * self.H0_cgs / (8. * np.pi * self.G_cgs)
        
        # critical density in units of [h2 Msun kpc-3]
        self.rho_crit0_cosmo = self.rho_crit0_cgs * (self.kpc3_cgs) / (self.Msun_cgs) / self.h_cosmo / self.h_cosmo

        # Normalization factors from Initialize Cosmology
        self.r0_gas = 1.0 # simulation ran with gas in [h-1 kpc] (???????)
        self.t0_gas = self.t_H0_cosmo / self.h_cosmo
        self.v0_gas = self.r0_gas / self.t0_gas
        self.rho0_gas = self.rho_crit0_cosmo * self.OmegaM
        self.phi0_gas = self.v0_gas * self.v0_gas  # energy units
        self.e0_gas = self.v0_gas * self.v0_gas
        self.p0_gas = self.rho0_gas * self.v0_gas * self.v0_gas # pressure units


        # conversion factors between cosmo [kpc/km Msun kyr] and cgs [cm gram sec] units
        # these factors DO NOT account for comoving units (ie, scale factor not incorporated)
        # multiplying array (in cgs units) by array_cgs2cosmo provides array in cosmo units
        # multiplying array (in cosmo units) by array_cosmo2cgs provides array in cgs units
        self.density_cgs2cosmo = self.rho_crit0_cosmo # [h2 Msun kpc-3]
        self.density_cosmo2cgs = 1. / self.density_cgs2cosmo

        self.velocity_cgs2cosmo = self.km_cgs # [km s-1]
        self.velocity_cosmo2cgs = 1. / self.velocity_cgs2cosmo

        self.mom_cgs2cosmo = self.density_cgs2cosmo * self.km_cgs # [h2 Msun kpc-3  km s-1]
        self.mom_cosmo2cgs = 1. / self.mom_cgs2cosmo

        self.energy_cgs2cosmo = self.density_cgs2cosmo * self.km_cgs * self.km_cgs # [h2 Msun kpc-3  km2 s-2]
        self.energy_cosmo2cgs = 1. / self.mom_cgs2cosmo


class ChollaSnapCosmologyHead:
    '''
    Cholla Snapshot Cosmology header object
        Serves as a header holding information that combines a ChollaCosmologyHead
            with a specific scale factor with the snapshot header object.
        
        Initialized with:
            snapHead (ChollaSnapHead): provides current redshift
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, snapHead, cosmoHead):
        self.a = snapHead.a
        self.cosmoHead = cosmoHead

        # calculate & attach current Hubble rate in [km s-1 Mpc-1] and [s-1]
        self.Hubble_cosmo = self.Hubble()
        self.Hubble_cgs = self.Hubble_cosmo * self.cosmoHead.km_cgs / self.cosmoHead.Mpc_cgs # in cgs [s-1]


    def Hubble(self):
        '''
        Return the current Hubble parameter

        Args:
            ...
        Returns:
            H (float): Hubble parameter (km/s/Mpc)
        '''

        a2 = self.a * self.a
        a3 = a2 * self.a
        a4 = a3 * self.a
        DE_factor = (self.a)**(-3. * (1. + self.cosmoHead.w0 + self.cosmoHead.wa))
        DE_factor *= np.exp(-3. * self.cosmoHead.wa * (1. - self.a))

        H0_factor = (self.cosmoHead.OmegaR / a4) + (self.cosmoHead.OmegaM / a3)
        H0_factor += (self.cosmoHead.OmegaK / a2) + (self.cosmoHead.OmegaL * DE_factor)

        return self.cosmoHead.H0 * np.sqrt(H0_factor)

    
    def dvHubble(self, dx):
        '''
        Return the Hubble flow through a cell

        Args:
            dx (float): comoving distance between cells (kpc)
        Returns:
            (float): Hubble flow over a cell (km/s)
        '''
        # convert [kpc] to [h-1 kpc]
        dx_h = dx / self.cosmoHead.h_cosmo

        dxh_cgs = dx_h * self.cosmoHead.kpc_cgs # h^-1 kpc * (#cm / kpc) =  h^-1 cm
        dxh_Mpc = dxh_cgs / self.cosmoHead.Mpc_cgs # h^-1 cm / (#cm / Mpc) = h^-1 Mpc
        
        # convert to physical length
        dxh_Mpc_phys = dxh_Mpc * self.a

        return self.Hubble() * dxh_Mpc_phys




class ChollaCosmoCalculator:
    '''
    Cholla Cosmological Calculator object
        Serves as a calculator for a cosmology at a specific scale factor.
        
        Initialized with:
            snapCosmoHead (ChollaSnapCosmologyHead): provides current redshift
            dims (tuple): size of data sets to act on
            dtype (np type): (optional) numpy precision to initialize output arrays 

    Values are returned in code units unless otherwise specified.
    '''

    def __init__(self, snapCosmoHead, dims, dtype=np.float32):
        self.snapCosmoHead = snapCosmoHead
        self.dims = dims
        self.dtype = dtype    

    def create_arr(self):
        '''
        Create and return an empty array
        
        Args:
            ...
        Returns:
            (arr): array of initialized dimensions and datatype
        '''

        return np.zeros(self.dims, dtype=self.dtype)

    def physical_length(self, length_comov):
        '''
        Calculate the physical length from a comoving length

        Args:
            length_comov (float): comoving length
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(length_comov.shape, self.dims)

        # initialize array with dims shape
        arr = self.create_arr()

        arr[:] = length_comov * self.snapCosmoHead.a

        return arr


    def physical_density(self, density_comov):
        '''
        Calculate the physical density from a comoving density

        Args:
            density_comov (float): comoving density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density_comov.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()

        a3 = self.snapCosmoHead.a * self.snapCosmoHead.a * self.snapCosmoHead.a
        arr[:] = density_comov / a3

        return arr


    def density_cosmo2cgs(self, density_cosmo):
        '''
        Convert the density saved in cosmological units of [h2 Msun kpc-3]
            to cgs units of [g cm-3]. With the large orders of magnitude
            involved, this calculation is completed in log-space

        Args:
            density_cosmo (float): density in cosmological units
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density_cosmo.shape, self.dims)

        # initialize array with dims shape
        arr = self.create_arr()
        
        # calculate h^2
        h_cosmo2 = self.snapCosmoHead.cosmoHead.h_cosmo * self.snapCosmoHead.cosmoHead.h_cosmo

        # take log of constants
        ln_hcosmo2 = np.log(h_cosmo2)
        ln_Msun = np.log(self.snapCosmoHead.cosmoHead.Msun_cgs)
        ln_kpc3 = np.log(self.snapCosmoHead.cosmoHead.kpc3_cgs)
        
        # take log of density
        ln_density_cosmo = np.log(density_cosmo)
        
        # convert values to cgs
        ln_density_cgs = ln_density_cosmo + ln_Msun + ln_hcosmo2 - ln_kpc3
        
        # take exp of log to get physical values
        arr[:] = np.exp(ln_density_cgs) # [g cm-3]

        return arr


    def velocity_cosmo2cgs(self, velocity_cosmo):
        '''
        Convert the velocity saved in cosmology units of [km s-1] to the cgs
            units of [cm s-1].

        Args:
            velocity_cosmo (float): velocity in cosmological units
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(velocity_cosmo.shape, self.dims)

        # initialize array with dims shape
        arr = self.create_arr()

        arr[:] = velocity_cosmo * self.snapCosmoHead.cosmoHead.km_cgs # [cm s-1]

        return arr


