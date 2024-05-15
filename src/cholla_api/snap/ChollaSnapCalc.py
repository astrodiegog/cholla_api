from cholla_api.snap.ChollaSnap import ChollaSnap

from cholla_api.data.ChollaHydroBox import ChollaHydroBox

from cholla_api.analysis.ChollaCalculator import ChollaCalculator


class ChollaSnapCalc:
    '''
    Cholla Snapshot Calculation object
        This object will handle the arguments that are passed onto 
            ChollaCalculator for some given SnapHead
        
        Will have two types of methods: get_calc and get_calcbox that allows
            user to create a calculation either on a specific box or on the
            entire sim global domain. get_calc methods will call 
            self.Snap.get_hydrodata which creates large arrays, so use with
            caution. get_calcbox will create a databox object and use its
            get_data method
    
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
        Calculate the hydro velocity magnitude for each requested nBox inside of
            nBoxes, and concatenate onto a global array for the entire sim.
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): array of the velocity magnitude
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        momx_str = self.Snap.head.DataHead.HydroHead.momx_str
        momy_str = self.Snap.head.DataHead.HydroHead.momy_str
        momz_str = self.Snap.head.DataHead.HydroHead.momz_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        # create global velocity magnitude array
        vmag_global = self.Calc.create_arr()
        
        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)
            
            density_box = box.get_data(namebase, dataDir, density_str, 
                                       old_format)
            momx_box = box.get_data(namebase, dataDir, momx_str, old_format)
            momy_box = box.get_data(namebase, dataDir, momy_str, old_format)
            momz_box = box.get_data(namebase, dataDir, momz_str, old_format)
            
            vmag_box = box_calc.velmag(density_box, momx_box, momy_box, momz_box)
            box.place_anydata(vmag_box, vmag_global)
            
        return vmag_global
    
    def get_vmag_box(self, namebase, dataDir, nBox):
        '''
        Calculate the hydro velocity magnitude for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            (arr): array of the velocity magnitude
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        momx_str = self.Snap.head.DataHead.HydroHead.momx_str
        momy_str = self.Snap.head.DataHead.HydroHead.momy_str
        momz_str = self.Snap.head.DataHead.HydroHead.momz_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        return box_calc.velmag(box.get_data(namebase, dataDir, density_str, 
                                            old_format),
                               box.get_data(namebase, dataDir, momx_str, 
                                            old_format),
                               box.get_data(namebase, dataDir, momy_str, 
                                            old_format),
                               box.get_data(namebase, dataDir, momz_str, 
                                            old_format))
    
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        # create global pressure array
        press_global = self.Calc.create_arr()
        
        # grab info needed to load in data
        old_format = self.Snap.head.DataHead.old_format
        
        if DE_flag:           
            # grab info needed to load in data
            gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
            
            for box_head in hydro_box_heads:
                # create a box + calculator object
                box = ChollaHydroBox(box_head, nsnap)
                box_calc = ChollaCalculator(box_head.local_dims)

                gasenergy_box = box.get_data(namebase, dataDir, gasenergy_str, 
                                             old_format)
                
                pressure_box = box_calc.pressure_DE(gasenergy_box, gamma)
                box.place_anydata(pressure_box, press_global)
                
        else:
            # grab info needed to load in data
            energy_str = self.Snap.head.DataHead.HydroHead.energy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            momx_str = self.Snap.head.DataHead.HydroHead.momx_str
            momy_str = self.Snap.head.DataHead.HydroHead.momy_str
            momz_str = self.Snap.head.DataHead.HydroHead.momz_str
            
            for box_head in hydro_box_heads:
                # create a box + calculator object
                box = ChollaHydroBox(box_head, nsnap)
                box_calc = ChollaCalculator(box_head.local_dims)

                energy_box = box.get_data(namebase, dataDir, energy_str, 
                                          old_format)
                density_box = box.get_data(namebase, dataDir, density_str, 
                                           old_format)
                momx_box = box.get_data(namebase, dataDir, momx_str, old_format)
                momy_box = box.get_data(namebase, dataDir, momy_str, old_format)
                momz_box = box.get_data(namebase, dataDir, momz_str, old_format)
                
                pressure_box = box_calc.pressure_noDE(energy_box, density_box, 
                                                      momx_box, momy_box, 
                                                      momz_box, gamma)
                
                box.place_anydata(pressure_box, press_global)
            
        return press_global
        
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        if DE_flag:
            # grab info needed to load in data
            gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
            
            return box_calc.pressure_DE(box.get_data(namebase, dataDir, 
                                                     gasenergy_str, old_format), 
                                        gamma)
        
        else:
            # grab info needed to load in data
            energy_str = self.Snap.head.DataHead.HydroHead.energy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            momx_str = self.Snap.head.DataHead.HydroHead.momx_str
            momy_str = self.Snap.head.DataHead.HydroHead.momy_str
            momz_str = self.Snap.head.DataHead.HydroHead.momz_str
            
            return box_calc.pressure_noDE(box.get_data(namebase, dataDir, 
                                                       energy_str, old_format),
                                          box.get_data(namebase, dataDir, 
                                                       density_str, old_format),
                                          box.get_data(namebase, dataDir, 
                                                       momx_str, old_format),
                                          box.get_data(namebase, dataDir, 
                                                       momy_str, old_format),
                                          box.get_data(namebase, dataDir, 
                                                       momz_str, old_format),
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        # create global internal energy arrays
        intenergy_global = self.Calc.create_arr()
        
        # grab info needed to load in data
        old_format = self.Snap.head.DataHead.old_format
        
        if DE_flag:           
            # grab info needed to load in data
            gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            
            for box_head in hydro_box_heads:
                # create a box + calculator object
                box = ChollaHydroBox(box_head, nsnap)
                box_calc = ChollaCalculator(box_head.local_dims)

                gasenergy_box = box.get_data(namebase, dataDir, gasenergy_str, 
                                             old_format)
                density_box = box.get_data(namebase, dataDir, density_str, 
                                           old_format)
                
                intenergy_box = box_calc.intenergy_DE(gasenergy_box, 
                                                      density_box, gamma)
                box.place_anydata(intenergy_box, intenergy_global)
                
        else:
            # grab info needed to load in data
            energy_str = self.Snap.head.DataHead.HydroHead.energy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            momx_str = self.Snap.head.DataHead.HydroHead.momx_str
            momy_str = self.Snap.head.DataHead.HydroHead.momy_str
            momz_str = self.Snap.head.DataHead.HydroHead.momz_str
            
            for box_head in hydro_box_heads:
                # create a box + calculator object
                box = ChollaHydroBox(box_head, nsnap)
                box_calc = ChollaCalculator(box_head.local_dims)

                energy_box = box.get_data(namebase, dataDir, energy_str, 
                                          old_format)
                density_box = box.get_data(namebase, dataDir, density_str, 
                                           old_format)
                momx_box = box.get_data(namebase, dataDir, momx_str, old_format)
                momy_box = box.get_data(namebase, dataDir, momy_str, old_format)
                momz_box = box.get_data(namebase, dataDir, momz_str, old_format)
                
                intenergy_box = box_calc.intenergy_noDE(energy_box, density_box, 
                                                       momx_box, momy_box, 
                                                       momz_box, gamma)
                
                box.place_anydata(intenergy_box, intenergy_global)
            
        return intenergy_global
    
    
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        if DE_flag:
            # grab info needed to load in data
            gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            
            return box_calc.intenergy_DE(box.get_data(namebase, dataDir, 
                                                      gasenergy_str, old_format), 
                                         box.get_data(namebase, dataDir, 
                                                      density_str, old_format),
                                         gamma)
        
        else:
            # grab info needed to load in data
            energy_str = self.Snap.head.DataHead.HydroHead.energy_str
            density_str = self.Snap.head.DataHead.HydroHead.density_str
            momx_str = self.Snap.head.DataHead.HydroHead.momx_str
            momy_str = self.Snap.head.DataHead.HydroHead.momy_str
            momz_str = self.Snap.head.DataHead.HydroHead.momz_str
            
            return box_calc.intenergy_noDE(box.get_data(namebase, dataDir, 
                                                        energy_str, old_format),
                                           box.get_data(namebase, dataDir, 
                                                        density_str, old_format),
                                           box.get_data(namebase, dataDir, 
                                                        momx_str, old_format),
                                           box.get_data(namebase, dataDir, 
                                                        momy_str, old_format),
                                           box.get_data(namebase, dataDir, 
                                                        momz_str, old_format),
                                           gamma)
    
    
    
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        # create global gas temperature arrays
        gastemp_global = self.Calc.create_arr()
        
        # grab info needed to load in data
        gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format

        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)

            gasenergy_box = box.get_data(namebase, dataDir, gasenergy_str, 
                                         old_format)
            density_box = box.get_data(namebase, dataDir, density_str, 
                                       old_format)

            gastemp_box = box_calc.gas_temp(gasenergy_box, density_box, gamma, 
                                            mu, velocity_unit)
            box.place_anydata(gastemp_box, gastemp_global)
            
        return gastemp_global
    
    
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
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        # grab info needed to load in data
        gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format

        return box_calc.gas_temp(box.get_data(namebase, dataDir, gasenergy_str, 
                                              old_format), 
                                 box.get_data(namebase, dataDir, density_str, 
                                              old_format),
                                 gamma, mu, velocity_unit)
    
    
        
    def get_overdensity(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the overdensity (density normalized by mean) for the entire
            sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D array of the overdensity
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        # create global overdensity arrays
        overdensity_global = self.Calc.create_arr()
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format

        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)

            density_box = box.get_data(namebase, dataDir, density_str, 
                                       old_format)

            overdensity_box = box_calc.overdensity_mean(density_box)
            box.place_anydata(overdensity_box, overdensity_global)
            
        return overdensity_global
        
        
    def get_overdensity_box(self, namebase, dataDir, nBox):
        '''
        Calculate the overdensity (density normalized by mean) for a subvolume 
            box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            (arr): 2D array of the overdensity
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format

        return box_calc.overdensity_mean(box.get_data(namebase, dataDir, density_str, 
                                                        old_format))
        
        
    def get_xprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the x-projected density for the entire sim box. Loop over 
            nBoxes, apply x projection, and place onto larger array using the 
            box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        new_dims = (self.Snap.head.DataHead.HydroHead.dims[1], 
                    self.Snap.head.DataHead.HydroHead.dims[2])
        
        # create global array
        xproj_global = self.Calc.create_subarr(new_dims)
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)
            
            # x projection is with k_index = 0
            xproj_box = box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                                  density_str, 
                                                                  old_format), 0)
            
            startX, startY, startZ = box_head.offset
            localX, localY, localZ = box_head.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            xproj_global[startY:endY, startZ:endZ] += xproj_box
        
        return xproj_global
    
    def get_xprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the x-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        # x projection is with k_index = 0
        return box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                         density_str, 
                                                         old_format), 0)
    
    def get_yprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the y-projected density for the entire sim box. Loop over 
            nBoxes, apply y projection, and place onto larger array using the 
            box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        new_dims = (self.Snap.head.DataHead.HydroHead.dims[0], 
                    self.Snap.head.DataHead.HydroHead.dims[2])
        
        # create global array
        yproj_global = self.Calc.create_subarr(new_dims)
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)
            
            # y projection is with k_index = 1
            yproj_box = box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                                  density_str, 
                                                                  old_format), 1)
            
            startX, startY, startZ = box_head.offset
            localX, localY, localZ = box_head.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            yproj_global[startX:endX, startZ:endZ] += yproj_box
        
        return yproj_global
    
    def get_yprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the y-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        # y projection is with k_index = 1
        return box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                         density_str, 
                                                         old_format), 1)
    
    
    def get_zprojection(self, namebase, dataDir, nBoxes=None):
        '''
        Calculate the z-projected density for the entire sim box. Loop over 
            nBoxes, apply y projection, and place onto larger array using the 
            box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        new_dims = (self.Snap.head.DataHead.HydroHead.dims[0], 
                    self.Snap.head.DataHead.HydroHead.dims[1])
        
        # create global array
        zproj_global = self.Calc.create_subarr(new_dims)
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab info needed to create ChollaHydroBoxes
        nsnap = self.Snap.head.nSnap
        
        # grab requested hydro box heads + ensure they have local_dims + offset
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        for box_head in hydro_box_heads:
            # create a box + calculator object
            box = ChollaHydroBox(box_head, nsnap)
            box_calc = ChollaCalculator(box_head.local_dims)
            
            # z projection is with k_index = 2
            zproj_box = box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                                  density_str, 
                                                                  old_format), 2)
            
            startX, startY, startZ = box_head.offset
            localX, localY, localZ = box_head.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            zproj_global[startX:endX, startY:endY] += zproj_box
        
        return yproj_global
    
    def get_zprojection_box(self, namebase, dataDir, nBox):
        '''
        Calculate the z-projected density for a subvolume box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.check_nbox(nBox)
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab info needed to load in data
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        old_format = self.Snap.head.DataHead.old_format
        
        # grab hydro box head
        hbox_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims before creating new calculator
        assert hbox_head.head_set
        box = ChollaHydroBox(hbox_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(hbox_head.local_dims)
        
        # z projection is with k_index = 2
        return box_calc.densityk_projection(box.get_data(namebase, dataDir, 
                                                         density_str, 
                                                         old_format), 2)

