from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.analysis.ChollaCalculator import ChollaCalculator


class ChollaSnapCalc:
    '''
    Cholla Snapshot Calculation object
        This object will handle the arguments that are passed onto 
            ChollaCalculator for some given SnapHead
        
        Will have two types of methods: get_calc and get_calcbox that allows
            user to create a calculation either on a specific box or on for the
            entire sim global domain
    
    TODO:
        once I have a method in ChollaSnap to load data within a subvolume,
            create a new method for each of these calculations
    '''
    
    def __init__(self, SnapHead):
        self.Snap = ChollaSnap(SnapHead)
        
        # use Hydro dims as global dims
        self.Calc = ChollaCalculator(self.Snap.head.DataHead.HydroHead.dims)
    
    def get_vmag(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the hydro velocity magnitude for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): array of the velocity magnitude
        '''
        
        return self.Calc.velmag(self.Snap.get_hydrodata(namebase, dataDir, 
                                                        "density", nBoxes),
                                self.Snap.get_hydrodata(namebase, dataDir, 
                                                        "momentum_x", nBoxes),
                                self.Snap.get_hydrodata(namebase, dataDir, 
                                                        "momentum_y", nBoxes),
                                self.Snap.get_hydrodata(namebase, dataDir, 
                                                        "momentum_z", nBoxes))
    
    def get_vmag_box(self, namebase, dataDir, nBox):
        '''
        Calculate the hydro velocity magnitude for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            nBox (int): index of the box to load
        Returns:
            (arr): array of the velocity magnitude
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.velmag(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                         "density", nBox),
                              self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                         "momentum_x", nBox),
                              self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                         "momentum_y", nBox),
                              self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                         "momentum_z", nBox))
    
    def get_pressure(self, namebase, dataDir, gamma, DE_flag, nBoxes=None):
        '''
        Calculate the pressure for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            DE_flag (bool): whether dual energy formalism was used
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): array of the pressure
        '''
        
        if DE_flag:
            return self.Calc.pressure_DE(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                 "GasEnergy", nBoxes), 
                                         gamma)
        else:
            return self.Calc.pressure_noDE(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "Energy", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "density", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_x", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_y", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_z", nBoxes),
                                           gamma)
        
    def get_pressure_box(self, namebase, dataDir, gamma, DE_flag, nBox):
        '''
        Calculate the pressure for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            DE_flag (bool): whether dual energy formalism was used
            nBox (int): index of the box to load
        Returns:
            (arr): array of the pressure
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        if DE_flag:
            return newCalc.pressure_DE(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                  "GasEnergy", nBox),
                                       gamma)
        else:
            return newCalc.pressure_noDE(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                    "Energy", nBox),
                                         self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                    "density", nBox),
                                         self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                    "momentum_x", nBox),
                                         self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                    "momentum_y", nBox),
                                         self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                    "momentum_z", nBox),
                                         gamma)
        
    def get_intenergy(self, namebase, dataDir, DE_flag, nBoxes=None):
        '''
        Calculate the specific internal energy for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            DE_flag (bool): whether dual energy formalism was used
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): array of the internal energy
        '''
        
        if DE_flag:
            return self.Calc.intenergy_DE(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                 "GasEnergy", nBoxes),
                                          self.Snap.get_hydrodata(namebase, dataDir, 
                                                                 "density", nBoxes))
        else:
            return self.Calc.intenergy_noDE(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "Energy", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "density", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_x", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_y", nBoxes),
                                           self.Snap.get_hydrodata(namebase, dataDir, 
                                                                   "momentum_z", nBoxes))

    def get_intenergy_box(self, namebase, dataDir, DE_flag, nBox):
        '''
        Calculate the specific internal energy for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            DE_flag (bool): whether dual energy formalism was used
            nBox (int): index of the box to load
        Returns:
            (arr): array of the internal energy
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        if DE_flag:
            return newCalc.intenergy_DE(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                   "GasEnergy", nBox),
                                        self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                   "density", nBox))
        else:
            return newCalc.intenergy_noDE(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "Energy", nBox),
                                          self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "density", nBox),
                                          self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "momentum_x", nBox),
                                          self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "momentum_y", nBox),
                                          self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "momentum_z", nBox))
        
    def get_gastemp(self, namebase, dataDir, gamma, mu, velocity_unit, nBoxes=None):
        '''
        Calculate the temperature of the gas for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): array of the temperature
        '''
        
        return self.Calc.gas_temp(self.Snap.get_hydrodata(namebase, dataDir, 
                                                          "GasEnergy", nBoxes),
                                  self.Snap.get_hydrodata(namebase, dataDir, 
                                                          "density", nBoxes), 
                                  gamma, mu, velocity_unit)
    
    def get_gastemp_box(self, namebase, dataDir, gamma, mu, velocity_unit, nBox):
        '''
        Calculate the temperature of the gas for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
            nBox (int): index of the box to load
        Returns:
            (arr): array of the temperature
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.gas_temp(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                           "GasEnergy", nBox),
                                self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                           "density", nBox), 
                                gamma, mu, velocity_unit)
        
    def get_overdensity(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the overdensity (density normalized by median) for the entire
            sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the overdensity
        '''
        
        return self.Calc.overdensity_median(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                    "density", nBoxes))
    
    def get_overdensity_box(self, namebase, dataDir, nBox):
        '''
        Calculate the overdensity (density normalized by median) for a subvolume 
            box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            (arr): 2D array of the overdensity
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.overdensity_median(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                     "density", nBox))
        
        
    def get_xprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the x-projected density for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        return self.Calc.densityk_projection(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                     "density", nBoxes),
                                             0)
    
    def get_xprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the x-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.densityk_projection(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                      "density", nBox),
                                             0)
    
    def get_yprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the y-projected density for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        return self.Calc.densityk_projection(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                     "density", nBoxes),
                                             1)
    
    def get_yprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the y-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.densityk_projection(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                      "density", nBox),
                                             1)
    
    def get_zprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the z-projected density for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        return self.Calc.densityk_projection(self.Snap.get_hydrodata(namebase, dataDir, 
                                                                     "density", nBoxes),
                                             2)
    
    def get_zprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the z-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the projected density
        '''
        
        HydroBoxHead = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert HydroBoxHead.head_set
        
        newCalc = ChollaCalculator(HydroBoxHead.local_dims)
        
        return newCalc.densityk_projection(self.Snap.get_hydroboxdata(namebase, dataDir, 
                                                                      "density", nBox),
                                           2)