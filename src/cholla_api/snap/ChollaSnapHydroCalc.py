import numpy as np

from cholla_api.analysis.ChollaHydroCalculator import ChollaHydroCalculator
from cholla_api.data.ChollaBox import ChollaBox
from cholla_api.data.ChollaBoxHydroCalc import ChollaBoxHydroCalc


class ChollaSnapHydroCalc:
    '''
    Cholla Snap Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - chSnap (ChollaSnap): ChollaSnap object
        - chGrid (ChollaGrid): ChollaGrid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
        - dtype (datatype): (optional) precision to calculate data
        - calc_box (bool): (optional) whether to use ChollaBoxHydroCalculator object

        if calc_box is set to True, then we use ChollaBox's place_data method
            to create a global data array
        else, then we instead complete calculation on the global primitive
            data. take care when setting calc_box=False because it can lead
            to loading large numpy arrays into memory !
    '''
    
    def __init__(self, chSnap, chGrid, chMacroFlags, dtype=np.float32, calc_box=True):
        self.Snap = ChollaSnap
        self.Grid = ChollaGrid
        self.Calculator = ChollaHydroCalculator((self.Grid.nx_global, self.Grid.ny_global, 
                                                 self.Grid.nz_global), 
                                                dtype=dtype)
        self.MacroFlags = chMacroFlags
        self.boxcalc = calc_box
        if self.boxcalc:
            from cholla_api.data.ChollaBoxHydroCalc import ChollaBoxHydroCalc
            # grab key strings from boxhead 0
            boxhead0 = self.Grid.get_BoxHead(0)
            box0 = ChollaBox(self.Snap.SnapPath, boxhead0)

    def get_vmag(self):
        '''
        Calculate and return the hydro velocity magnitude

        Args:
            ...
        Returns:
            (arr): velocity magnitude
        '''

        if self.boxcalc:
            vmag_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                vmag_box = boxhydrocalc.get_vmag()
                box.place_data(vmag_box, vmag_snap)

            return vmag_snap
        else:
            density_str = self.box0.density_str
            momx_str = self.box0.momx_str
            momy_str = self.box0.momy_str
            momz_str = self.box0.momz_str

            return self.Calculator.velmag(self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                  density_str, self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                  momx_str, self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                  momy_str, self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                  momz_str, self.Calculator.dtype))

    def get_intenergy(self):
        '''
        Calculate and return the internal energy

        Args:
            ...
        Returns:
            (arr): internal energy
        '''

        if self.boxcalc:
            intenergy_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags,
                                                  self.Calculator.dtype)
                intenergy_box = boxhydrocalc.get_intenergy()
                box.place_data(intenergy_box, vmag_snap)

            return intenergy_snap
    
        else:
            density_str = self.Box.density_str
            energy_str = self.Box.energy_str
            momx_str = self.Box.momx_str
            momy_str = self.Box.momy_str
            momz_str = self.Box.momz_str

            return self.Calculator.int_energy(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                      energy_str, self.Calculator.dtype),
                                              self.Snap.get_hydrodata(self.Grid, self.MacroFlags
                                                                      density_str, self.Calculator.dtype),
                                              self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                      momx_str, self.Calculator.dtype),
                                              self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                      momy_str, self.Calculator.dtype),
                                              self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                      momz_str, self.Calculator.dtype))


    def get_pressure(self, gamma):
        '''
        Calculate and return the pressure

        Args:
            gamma (float): ratio of specific heats
        Returns:
            (arr): pressure
        '''

        if self.boxcalc:
            pressure_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                pressure_box = boxhydrocalc.get_pressure(gamma)
                box.place_data(pressure_box, pressure_snap)

            return pressure_snap
        else:
            if self.MacroFlags.DualEnergy:
                gasenergy_str = self.box0.gasenergy_str

                return self.Calculator.pressure_DE(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                           gasenergy_str, self.Calculator.dtype),
                                                   gamma)
            else:
                energy_str = self.box0.energy_str
                density_str = self.box0.density_str
                momx_str = self.box0.momx_str
                momy_str = self.box0.momy_str
                momz_str = self.box0.momz_str

                return self.Calculator.pressure_noDE(self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                             energy_str, self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                             density_str, self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                             momx_str, self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                             momy_str, self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                             momz_str, self.Calculator.dtype),
                                                     gamma)


    def get_specintenergy(self):
        '''
        Calculate and return the specific internal energy

        Args:
            ...
        Returns:
            (arr): internal energy
        '''

        if self.boxcalc:
            specintenergy_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                specintenergy_box = boxhydrocalc.get_specintenergy()
                box.place_data(specintenergy_box, specintenergy_snap)

            return pressure_snap
        else:
            density_str = self.box0.density_str
            if self.MacroFlags.DualEnergy:
                gasenergy_str = self.box0.gasenergy_str

                return self.Calculator.specintenergy_DE(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                                gasenergy_str, self.Calculator.dtype),
                                                        self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                                density_str, self.Calculator.dtype)) 
            else:
                energy_str = self.box0.energy_str
                momx_str = self.box0.momx_str
                momy_str = self.box0.momy_str
                momz_str = self.box0.momz_str

                return self.Calculator.specintenergy_noDE(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                                  energy_str, self.Calculator.dtype),
                                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                                  density_str, self.Calculator.dtype),
                                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                                  momx_str, self.Calculator.dtype),
                                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                                  momy_str, self.Calculator.dtype),
                                                          self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                                  momz_str, self.Calculator.dtype))

    def get_gastemp(self, gamma, mu, energy_unit):
        '''
        Calculate and return the gas temperature

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            energy_unit (float): unit conversion from energy code units to cgs
        Returns:
            (arr): gas temperature
        '''

        if self.boxcalc:
            gastemp_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                gastemp_box = boxhydrocalc.get_gastemp(gamma, mu, energy_unit)
                box.place_data(gastemp_box, gastemp_snap)

            return gastemp_snap
        else:
            gasenergy_str = self.box0.gasenergy_str
            density_str = self.box0.density_str

            return self.Calculator.gastemp(self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                   gasenergy_str, self.Calculator.dtype),
                                           self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                   density_str, self.Calculator.dtype),
                                           gamma, mu, energy_unit)
        
    def get_overdensity(self):
        '''
        Calculate the density normalized by the mean density

        Args:
            ...
        Returns:
            (arr): overdensity
        '''

        if self.boxcalc:
            overdensity_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                overdensity_box = boxhydrocalc.get_overdensity()
                box.place_data(overdensity_box, overdensity_snap)

            return overdensity_snap
        else:
            density_str = self.box0.density_str

            return self.Calculator.overdensity_mean(self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                            density_str, self.Calculator.dtype))
       
    def get_xprojection(self):
        '''
        Calculate the x-projected density

        Args:
            ...
        Returns:
            (arr): x-projected density
        '''
        
        if self.boxcalc:
            nx_snap = self.Calculator.create_subarr((self.Grid.ny_global, 
                                                     self.Grid.nz_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                nx_box = boxhydrocalc.get_xprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                nx_snap[startY:endY, startZ:endZ] += nx_box

            return nx_snap
        else:
            density_str = self.box0.density_str

            # x-projection uses k-index=0
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                               density_str, self.Calculator.dtype), 
                                                       0)
    
    def get_yprojection(self):
        '''
        Calculate the y-projected density

        Args:
            ...
        Returns:
            (arr): y-projected density
        '''

        if self.boxcalc:
            ny_snap = self.Calculator.create_subarr((self.Grid.nx_global,
                                                     self.Grid.nz_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags, 
                                                  self.Calculator.dtype)
                ny_box = boxhydrocalc.get_yprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                ny_snap[startX:endX, startZ:endZ] += ny_box

            return ny_snap
        else:
            density_str = self.box0.density_str

            # y-projection uses k-index=1
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, self.MacroFlags, 
                                                                               density_str, self.Calculator.dtype),
                                                       1)

    def get_zprojection(self):
        '''
        Calculate the z-projected density

        Args:
            ...
        Returns:
            (arr): z-projected density
        '''

        if self.boxcalc:
            nz_snap = self.Calculator.create_subarr((self.Grid.nx_global,
                                                     self.Grid.ny_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags,
                                                  self.Calculator.dtype)
                nz_box = boxhydrocalc.get_zprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                nz_snap[startX:endX, startY:endY] += nz_box

            return nz_snap
        else:
            density_str = self.box0.density_str

            # z-projection uses k-index=2
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, self.MacroFlags,
                                                                               density_str, self.Calculator.dtype),
                                                       2)



    def get_phasespace(self, gamma, mu, energy_unit):
        '''
        Calculate and return the phase space

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            energy_unit (float): unit conversion from energy code units to cgs
        Returns:
            (tuple): phase space, overdensity bin edges, temp bin edges
        '''

        if self.boxcalc:
            phasespace_snap = self.Calculator.create_subarr((49,49))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead, self.MacroFlags)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.MacroFlags,
                                                  self.Calculator.dtype)
                phasespace_box, xbin, ybin  = boxhydrocalc.get_phasespace(gamma, mu, energy_unit)
                phasespace_snap += phasespace_box

            return (phasespace_snap, xbin, ybin)
        else:
            return self.Calculator.create_phase(np.log10(self.get_gastemp(gamma, mu, 
                                                                          energy_unit)).ravel(),
                                                np.log10(self.get_overdensity()).ravel())


