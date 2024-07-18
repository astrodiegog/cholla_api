
import numpy as np

from cholla_api.analysis.ChollaHydroCalculator import ChollaHydroCalculator
from cholla_api.data.ChollaBox import ChollaBox


class ChollaSnapHydroCalc:
    '''
    Cholla Snap Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - ChollaSnap (ChollaSnap): ChollaSnap object
        - ChollaGrid (ChollaGrid): ChollaGrid object
        - dtype (datatype): precision to calculate data
        - calc_box (bool): whether to use ChollaBoxHydroCalculator object

        if calc_box is set to True, then we use ChollaBox's place_data method
            to create a global data array
        else, then we instead complete calculation on the global primitive
            data. take care when setting calc_box=False because it can lead
            to loading large numpy arrays into memory !
    '''
    
    def __init__(self, ChollaSnap, ChollaGrid, dtype=np.float32, calc_box=True):
        self.Snap = ChollaSnap
        self.Grid = ChollaGrid
        self.Calculator = ChollaHydroCalculator((self.Grid.nx_global, self.Grid.ny_global, 
                                                 self.Grid.nz_global), 
                                                dtype=dtype)
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
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                vmag_box = boxhydrocalc.get_vmag()
                box.place_data(vmag_box, vmag_snap)

            return vmag_snap
        else:
            density_str = self.box0.density_str
            momx_str = self.box0.momx_str
            momy_str = self.box0.momy_str
            momz_str = self.box0.momz_str

            return self.Calculator.velmag(self.Snap.get_hydrodata(self.Grid, density_str, 
                                                                  self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, momx_str, 
                                                                  self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, momy_str, 
                                                                  self.Calculator.dtype),
                                          self.Snap.get_hydrodata(self.Grid, momz_str, 
                                                                  self.Calculator.dtype))


    def get_pressure(self, DE_flag, gamma):
        '''
        Calculate and return the pressure

        Args:
            DE_flag (bool): whether to use dual-energy formalism
            gamma (float): ratio of specific heats
        Returns:
            (arr): pressure
        '''

        if self.boxcalc:
            pressure_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                pressure_box = boxhydrocalc.get_pressure(DE_flag, gamma)
                box.place_data(pressure_box, pressure_snap)

            return pressure_snap
        else:
            if DE_flag:
                gasenergy_str = self.box0.gasenergy_str

                return self.Calculator.pressure_DE(self.Snap.get_hydrodata(self.Grid, gasenergy_str,
                                                                           self.Calculator.dtype),
                                                   gamma)
            else:
                energy_str = self.box0.energy_str
                density_str = self.box0.density_str
                momx_str = self.box0.momx_str
                momy_str = self.box0.momy_str
                momz_str = self.box0.momz_str

                return self.Calculator.pressure_noDE(self.Snap.get_hydrodata(self.Grid, energy_str,
                                                                             self.Calculator.dtype),
                                                    self.Snap.get_hydrodata(self.Grid, density_str,
                                                                             self.Calculator.dtype),
                                                    self.Snap.get_hydrodata(self.Grid, momx_str,
                                                                             self.Calculator.dtype),
                                                    self.Snap.get_hydrodata(self.Grid, momy_str,
                                                                             self.Calculator.dtype),
                                                    self.Snap.get_hydrodata(self.Grid, momz_str,
                                                                             self.Calculator.dtype),
                                                    gamma)


    def get_intenergy(self, DE_flag):
        '''
        Calculate and return the internal energy

        Args:
            DE_flag (bool): whether to use dual-energy formalism
        Returns:
            (arr): internal energy
        '''

        if self.boxcalc:
            intenergy_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                intenergy_box = boxhydrocalc.get_intenergy(DE_flag, gamma)
                box.place_data(intenergy_box, intenergy_snap)

            return pressure_snap
        else:
            density_str = self.box0.density_str
            if DE_flag:
                gasenergy_str = self.box0.gasenergy_str

                return self.Calculator.intenergy_DE(self.Snap.get_hydrodata(self.Grid, gasenergy_str,
                                                                            self.Calculator.dtype),
                                                    self.Snap.get_hydrodata(self.Grid, density_str,
                                                                            self.Calculator.dtype)) 
            else:
                energy_str = self.box0.energy_str
                momx_str = self.box0.momx_str
                momy_str = self.box0.momy_str
                momz_str = self.box0.momz_str

                return self.Calculator.intenergy_noDE(self.Snap.get_hydrodata(self.Grid, energy_str,
                                                                              self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, density_str,
                                                                              self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, momx_str,
                                                                              self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, momy_str,
                                                                              self.Calculator.dtype),
                                                     self.Snap.get_hydrodata(self.Grid, momz_str,
                                                                              self.Calculator.dtype))

    def get_gastemp(self, gamma, mu, velocity_unit):
        '''
        Calculate and return the gas temperature

        Args:
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
        Returns:
            (arr): gas temperature
        '''

        if self.boxcalc:
            gastemp_snap = self.Calculator.create_arr()

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                gastemp_box = boxhydrocalc.get_gastemp(gamma, mu, velocity_unit)
                box.place_data(gastemp_box, gastemp_snap)

            return gastemp_snap
        else:
            gasenergy_str = self.box0.gasenergy_str
            density_str = self.box0.density_str

            return self.Calculator.gastemp(self.Snap.get_hydrodata(self.Grid, gasenergy_str,
                                                                   self.Calculator.dtype),
                                           self.Snap.get_hydrodata(self.Grid, density_str,
                                                                   self.Calculator.dtype),
                                           gamma, mu, velocity_unit)
        
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
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                overdensity_box = boxhydrocalc.get_overdensity()
                box.place_data(overdensity_box, overdensity_snap)

            return overdensity_snap
        else:
            density_str = self.box0.density_str

            return self.Calculator.overdensity_mean(self.Snap.get_hydrodata(self.Grid, density_str,
                                                                            self.Calculator.dtype))
       
    def get_xprojection(self):
        '''
        Calculate the x-projected density for Box

        Args:
            ...
        Returns:
            (arr): x-projected density
        '''
        
        if self.boxcalc:
            nx_snap = self.Calculator.create_arr((self.Grid.ny_global, 
                                                  self.Grid.nz_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                nx_box = boxhydrocalc.get_xprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                nx_snap[startY:endY, startZ:endZ] += nx_box

            return nx_snap
        else:
            density_str = self.box0.density_str

            # x-projection uses k-index=0
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, density_str,
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

        if self.boxcalc:
            ny_snap = self.Calculator.create_arr((self.Grid.nx_global,
                                                  self.Grid.nz_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                ny_box = boxhydrocalc.get_yprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                ny_snap[startX:endX, startZ:endZ] += ny_box

            return ny_snap
        else:
            density_str = self.box0.density_str

            # y-projection uses k-index=1
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, density_str,
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

        if self.boxcalc:
            nz_snap = self.Calculator.create_arr((self.Grid.nx_global,
                                                  self.Grid.ny_global))

            for boxhead in self.Grid.get_BoxHeads():
                box = ChollaBox(self.Snap.SnapPath, boxhead)
                boxhydrocalc = ChollaBoxHydroCalc(box, self.Calculator.dtype)
                nz_box = boxhydrocalc.get_zprojection()

                startX, startY, startZ = box.BoxHead.offset
                localX, localY, localZ = box.BoxHead.local_dims
                endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
                nz_snap[startX:endX, startY:endY] += nz_box

            return nz_snap
        else:
            density_str = self.box0.density_str

            # z-projection uses k-index=2
            return self.Calculator.densityk_projection(self.Snap.get_hydrodata(self.Grid, density_str,
                                                                               self.Calculator.dtype),
                                                       2)

    


