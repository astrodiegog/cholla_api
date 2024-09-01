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


class ChollaCosmoCalculator:
    '''
    Cholla Cosmological Calculator object
        Serves as a calculator where the calculated values have some expected
            size and datatype (default is float). Assert that inputs are of same
            shape as dims that was used to initialize this calculator. To
            complete some analysis, this ChollaCalculator will be the mediator
            that will act on the primitive saved values. Acts on a specific
            snapshot to convert between physical and comoving units

        Initialized with:
            snapHead (ChollaSnapHead): provides current redshift
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units

    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, snapHead, cosmoHead):
        self.snapHead = snapHead
        self.cosmoHead = cosmoHead

    def Hubble(self):
        '''
        Return the current Hubble parameter

        Args:
            ...
        Returns:
            H (float): Hubble parameter (km/s/Mpc)
        '''

        a2 = self.snapHead.a * self.snapHead.a
        a3 = a2 * self.snapHead.a
        a4 = a3 * self.snapHead.a
        DE_factor = (self.snapHead.a)**(-3. * (1. + self.cosmoHead.w0 + self.cosmoHead.wa))
        DE_factor *= np.exp(-3. * self.cosmoHead.wa * (1. - self.snapHead.a))

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
        dxh_Mpc_phys = self.physical_length(dxh_Mpc)

        return self.Hubble() * dxh_Mpc_phys


    def physical_length(self, length_comov):
        '''
        Calculate the physical length from a comoving length

        Args:
            length_comov (float): comoving length
        Returns:
            length_phys (float): physical length
            arr (arr): array that will hold data
        '''

        length_phys = length_comov * self.snapHead.a

        return length_comov


    def physical_density(self, density_comov):
        '''
        Calculate the physical density from a comoving density

        Args:
            density_comov (float): comoving density
        Returns:
            density_phys (float): physical density
        '''

        a3 = self.snapHead.a * self.snapHead.a * self.snapHead.a
        density_phys = density_comov / a3

        return density_comov


