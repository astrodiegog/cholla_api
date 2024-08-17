
import numpy as np

from cholla_api.analysis.ChollaHydroCalculator import ChollaHydroCalculator

class ChollaBoxHydroCalc:
    '''
    Cholla Box Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - ChollaBox (ChollaBox): ChollaBox object
        - dtype (datatype): precision to compute data
    '''
    
    def __init__(self, ChollaBox, dtype=np.float32):
        self.Box = ChollaBox
        self.Calculator = ChollaHydroCalculator(self.Box.BoxHead.local_dims, dtype=dtype)

    def get_vmag(self):
        '''
        Calculate and return the hydro velocity magnitude for cells in Box

        Args:
            ...
        Returns:
            (arr): velocity magnitude
        '''

        density_str = self.Box.density_str
        momx_str = self.Box.momx_str
        momy_str = self.Box.momy_str
        momz_str = self.Box.momz_str

        return self.Calculator.velmag(self.Box.get_hydrodata(density_str, 
                                                            self.Calculator.dtype),
                                      self.Box.get_hydrodata(momx_str,
                                                            self.Calculator.dtype),
                                      self.Box.get_hydrodata(momy_str,
                                                            self.Calculator.dtype),
                                      self.Box.get_hydrodata(momz_str,
                                                            self.Calculator.dtype))

    def get_pressure(self, DE_flag, gamma):
        '''
        Calculate and return the pressure for cells in Box

        Args:
            DE_flag (bool): whether to use dual-energy formalism
            gamma (float): ratio of specific heats
        Returns:
            (arr): pressure
        '''

        if DE_flag:
            gasenergy_str = self.Box.gasenergy_str

            return self.Calculator.pressure_DE(self.Box.get_hydrodata(gasenergy_str,
                                                                      self.Calculator.dtype),
                                               gamma)
        else:
            density_str = self.Box.density_str
            energy_str = self.Box.energy_str
            momx_str = self.Box.momx_str
            momy_str = self.Box.momy_str
            momz_str = self.Box.momz_str

            return self.Calculator.pressure_noDE(self.Box.get_hydrodata(energy_str,
                                                                        self.Calculator.dtype),
                                                 self.Box.get_hydrodata(density_str,
                                                                        self.Calculator.dtype),
                                                 self.Box.get_hydrodata(momx_str,
                                                                        self.Calculator.dtype),
                                                 self.Box.get_hydrodata(momy_str,
                                                                        self.Calculator.dtype),
                                                 self.Box.get_hydrodata(momz_str,
                                                                        self.Calculator.dtype),
                                                 gamma)


    def get_intenergy(self, DE_flag):
        '''
        Calculate and return the internal energy for cells in Box

        Args:
            DE_flag (bool): whether to use dual-energy formalism
        Returns:
            (arr): internal energy
        '''

        if DE_flag:
            gasenergy_str = self.Box.gasenergy_str
            density_str = self.Box.density_str

            return self.Calculator.intenergy_DE(self.Box.get_hydrodata(gasenergy_str,
                                                                       self.Calculator.dtype),
                                                self.Box.get_hydrodata(density_str,
                                                                       self.Calculator.dtype))
        else:
            density_str = self.Box.density_str
            energy_str = self.Box.energy_str
            momx_str = self.Box.momx_str
            momy_str = self.Box.momy_str
            momz_str = self.Box.momz_str

            return self.Calculator.intenergy_noDE(self.Box.get_hydrodata(energy_str,
                                                                         self.Calculator.dtype),
                                                  self.Box.get_hydrodata(density_str,
                                                                         self.Calculator.dtype),
                                                  self.Box.get_hydrodata(momx_str,
                                                                         self.Calculator.dtype),
                                                  self.Box.get_hydrodata(momy_str,
                                                                         self.Calculator.dtype),
                                                  self.Box.get_hydrodata(momz_str,
                                                                         self.Calculator.dtype))

    def get_gastemp(self, gamma, mu, energy_unit):
        '''
        Calculate and return the gas temperature for cells in Box

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            energy_unit (float): unit conversion from energy code units to cgs
        Returns:
            (arr): gas temperature
        '''
        
        gasenergy_str = self.Box.gasenergy_str
        density_str = self.Box.density_str

        return self.Calculator.gas_temp(self.Box.get_hydrodata(gasenergy_str,
                                                               self.Calculator.dtype),
                                        self.Box.get_hydrodata(density_str,
                                                               self.Calculator.dtype),
                                        gamma, mu, energy_unit)

    def get_overdensity(self):
        '''
        Calculate the density normalized by the mean density in Box

        Args:
            ...
        Returns:
            (arr): overdensity
        '''

        density_str = self.Box.density_str

        return self.Calculator.overdensity_mean(self.Box.get_hydrodata(density_str,
                                                                       self.Calculator.dtype))

    def get_xprojection(self):
        '''
        Calculate the x-projected density for Box

        Args:
            ...
        Returns:
            (arr): x-projected density
        '''
        
        density_str = self.Box.density_str

        # x-projection uses k-index=0
        return self.Calculator.densityk_projection(self.Box.get_hydrodata(density_str,
                                                                          self.Calculator.dtype),
                                                   0)
    
    def get_yprojection(self):
        '''
        Calculate the y-projected density for Box

        Args:
            ...
        Returns:
            (arr): y-projected density
        '''

        density_str = self.Box.density_str

        # y-projection uses k-index=1
        return self.Calculator.densityk_projection(self.Box.get_hydrodata(density_str,
                                                                          self.Calculator.dtype),
                                                   1)

    def get_zprojection(self):
        '''
        Calculate the z-projected density for Box

        Args:
            ...
        Returns:
            (arr): z-projected density
        '''

        density_str = self.Box.density_str

        # z-projection uses k-index=2
        return self.Calculator.densityk_projection(self.Box.get_hydrodata(density_str,
                                                                          self.Calculator.dtype),
                                                   2)
    
    def get_phasespace(self, gamma, mu, energy_unit):
        '''
        Calculate and return phase space for Box

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            energy_unit (float): unit conversion from energy code units to cgs
        Returns:
            (tuple): phase space, overdensity bin edges, temp bin edges
        '''

        return self.Calculator.create_phase(np.log10(self.get_gastemp(gamma, mu, 
                                                                      energy_unit)).ravel() ,
                                            np.log10(self.get_overdensity()).ravel() )



