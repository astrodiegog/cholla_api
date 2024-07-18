import numpy as np

class ChollaHydroCalculator:
    '''
    Cholla Calculator object
        Serves as a calculator where the calculated values have some expected
            size and datatype (default is float). Assert that inputs are of same
            shape as dims that was used to initialize this calculator. To 
            complete some analysis, this ChollaCalculator will be the mediator 
            that will act on the primitive saved values. 
        
        Contains methods to calculate
            hydro velx, vely, velz, velmag
            pressure
            specific internal energy
            temperature
    
    Values are returned in code units unless otherwise specified.
    '''
    def __init__(self, dims, dtype=np.float32):
        self.dims = dims
        self.dtype = dtype
        
        # cgs constant values
        self.mp = 1.672622e-24
        self.kB = 1.380658e-16
    
    def create_arr(self):
        '''
        Create and return an empty array
        
        Args:
            ...
        Returns:
            (arr): array of initialized dimensions and datatype
        '''
        
        return np.zeros(self.dims, dtype=self.dtype)
    
    def create_subarr(self, new_dims):
        '''
        Create and return an empty array of new dimension shapes
        
        Args:
            new_dims (tuple): shape of new array
        Returns:
            (arr): array of new dimensions and initialized datatype
        '''
        
        return np.zeros(new_dims, dtype=self.dtype)
    
    def velx(self, density, mom_x):
        '''
        Calculate the x-projected hydro velocity data
        
        Args:
            density (arr): hydrodynamic mass density
            mom_x (arr): x-momentum density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_x.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = mom_x/density
        
        return arr
    
    def vely(self, density, mom_y):
        '''
        Calculate the y-projected hydro velocity data
        
        Args:
            density (arr): hydrodynamic mass density
            mom_y (arr): y-momentum density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_y.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = mom_y/density
        
        return arr
    
    def velz(self, density, mom_z):
        '''
        Calculate the z-projected hydro velocity data
        
        Args:
            density (arr): hydrodynamic mass density
            mom_z (arr): z-momentum density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_z.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = mom_z/density
        
        return arr
    
    def velmag(self, density, mom_x, mom_y, mom_z):
        '''
        Calculate the hydro velocity magnitude
        
        Args:
            density (arr): hydrodynamic mass density
            mom_x (arr): x-momentum density
            mom_y (arr): y-momentum density
            mom_z (arr): z-momentum density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_x.shape, self.dims)
        assert np.array_equal(mom_y.shape, self.dims)
        assert np.array_equal(mom_z.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = np.sqrt((mom_x/density)**2 + (mom_y/density)**2 + (mom_z/density)**2)
        # arr[:] = np.sqrt(mom_x**2 + mom_y**2 + mom_z**2)/(density)
        
        return arr
    
    
    def pressure_DE(self, gas_energy, gamma):
        '''
        Calculate the pressure data for the dual-energy formalism
        
        Args:
            gas_energy (arr): the thermal energy density, equivalent to the 
                total energy density minus the kinetic energy density
            gamma (float): ratio of specific heats
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(gas_energy.shape, self.dims)
        assert gamma > 1.0
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = gas_energy*(gamma - 1.0)
        
        return arr
    
    def pressure_noDE(self, energy, density, mom_x, mom_y, mom_z, gamma):
        '''
        Calculate the pressure data without the dual-energy formalism
        
        Args:
            energy (arr): the total energy density
            density (arr): hydrodynamic mass density
            mom_x (arr): x-momentum density
            mom_y (arr): y-momentum density
            mom_z (arr): z-momentum density
            gamma (float): ratio of specific heats
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(energy.shape, self.dims)
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_x.shape, self.dims)
        assert np.array_equal(mom_y.shape, self.dims)
        assert np.array_equal(mom_z.shape, self.dims)
        assert gamma > 1.0
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        # could calculate vmag inline instead of calling function
        arr[:] = (energy - (0.5)*(density)*(self.velmag(density, mom_x, mom_y, mom_z))**2)*(gamma - 1.0)
        
        return arr
    
    
    def intenergy_DE(self, gas_energy, density):
        '''
        Calculate the specific internal energy data with the dual-energy 
            formalism
        
        Args:
            gas_energy (arr): the thermal energy density, equivalent to the 
                total energy density minus the kinetic energy density
            density (arr): hydrodynamic mass density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(gas_energy.shape, self.dims)
        assert np.array_equal(density.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = gas_energy/density
        
        return arr
    
    
    def intenergy_noDE(self, energy, density, mom_x, mom_y, mom_z):
        '''
        Calculate the specific internal energy data without the dual-energy 
            formalism
        
        Args:
            energy (arr): the total energy density
            density (arr): hydrodynamic mass density
            mom_x (arr): x-momentum density
            mom_y (arr): y-momentum density
            mom_z (arr): z-momentum density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(energy.shape, self.dims)
        assert np.array_equal(density.shape, self.dims)
        assert np.array_equal(mom_x.shape, self.dims)
        assert np.array_equal(mom_y.shape, self.dims)
        assert np.array_equal(mom_z.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        # like with pressure, vmag can be calculated in-line
        arr[:] = (energy/density) - (0.5)*(self.velmag(density, mom_x, mom_y, mom_z))**2
        
        return arr
    
    def gas_temp(self, gas_energy, density, gamma, mu, velocity_unit):
        '''
        Calculate the gas temperature in Kelvin
            We use the ideal gas equation kBT = (gamma - 1) * epsilon * mp * mu
            where epsilon is the specific internal energy, mp is the mass of 
            a proton, and mu is the mean molecular weight in amu.
            
            We note that epsilon is in units of energy per mass which is the 
            same as velocity units squared (for KE, mv^2 divided by mass is v^2)
            
            In effect, we need a conversion factor from code units to cm^2/s^2
        
        Args:
            gas_energy (arr): the thermal energy density, equivalent to the 
                total energy density minus the kinetic energy density
            density (arr): hydrodynamic mass density
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(gas_energy.shape, self.dims)
        assert np.array_equal(density.shape, self.dims)
        assert gamma > 1.0
        assert mu > 0.0
        assert velocity_unit > 0.0
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        coeff = (velocity_unit**2) * (gamma - 1.0) * self.mp * mu / self.kB
        arr[:] += gas_energy
        arr[:] *= coeff
        arr[:] /= density
        
        return arr
    
    
    def overdensity_median(self, density):
        '''
        Calculate the density normalized by the median of the density
        
        Args:
            density (arr): hydrodynamic mass density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = density/np.median(density)
        
        return arr
    
    
    def overdensity_mean(self, density):
        '''
        Calculate the density normalized by the mean of the density
        
        Args:
            density (arr): hydrodynamic mass density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        
        # initialize array with dims shape
        arr = self.create_arr()
        
        arr[:] = density/np.mean(density)
        
        return arr
    
    
    def densityk_projection(self, density, k_index):
        '''
        Calculate the density projection along some k_index.
            x-projection (y,z): k_index = 0
            y-projection (x,z): k_index = 1
            z-projection (x,y): k_index = 2
        
        Args:
            density (arr): hydrodynamic mass density
        Returns:
            arr (arr): array that will hold data
        '''
        assert np.array_equal(density.shape, self.dims)
        assert k_index >= 0
        assert k_index < 3
        
        # new shape will be of the other two indices, we don't assume cube dims
        # use a boolean mask of the shape to grab the other two dim sizes
        newshape_mask = np.ones(3, dtype=bool)
        newshape_mask[k_index] = False
        # need to cast density.shape tuple to an array to apply boolean mask
        subarr_shape = np.array(density.shape)[newshape_mask]
        
        # initialize array with dims shape
        arr = self.create_subarr(subarr_shape)
        
        arr[:] = np.sum(density, axis=k_index)
        
        return arr
 

    def create_phase(self, log_temp, log_overdensity):
        '''
        Create a 2D histogram for temperature and overdensity. Designed for
            cosmological simulations
        
        Args:
            log_temp (arr): 1D log-base10 of the gas temperature
            log_temp (arr): 1D log-base10 of the gas overdensity
        Returns:
            (tuple): 2D histogram, bin edges along 1st dimension, and bin edges 
                along 2nd dimension
        '''

        assert log_temp.ndim == 1
        assert log_overdensity.ndim == 1

        min_loverdensity, max_loverdensity = -2, 4
        min_ltemp, max_ltemp = 2, 8
        # min_ltemp, max_ltemp = 8, 14

        phasespace = np.histogram2d(log_overdensity, log_temp,
                                    range=((min_loverdensity, max_loverdensity),
                                           (min_ltemp, max_ltemp)),
                                    bins=(49,49))

        return (phasespace, log_overdensity_bins, log_temp_bins)

