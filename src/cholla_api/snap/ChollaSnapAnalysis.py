import numpy as np
from time import time

from cholla_api.data.ChollaHydroBox import ChollaHydroBox
from cholla_api.snap.ChollaSnap import ChollaSnap


from cholla_api.analysis.ChollaCalculator import ChollaCalculator
from cholla_api.analysis.ChollaAnalysis import ChollaAnalysis
from cholla_api.analysis.ChollaVizAnalysis import ChollaVizAnalysis

class ChollaSnapAnalysis:
    '''
    Cholla Snapshot Analysis object
        This object will handle the arguments that are passed onto 
            ChollaAnalysis for some given SnapHead
        
        Similar to SnapCalc, where we will have global and nBox specific methods
    '''
    
    def __init__(self, SnapHead):
        self.Snap = ChollaSnap(SnapHead)
        self.Calc = ChollaCalculator(self.Snap.head.DataHead.HydroHead.dims)
        self.Analysis = ChollaAnalysis(self.Snap.head.DataHead.HydroHead.dims)
    
    
    def get_phasespace_box(self, namebase, dataDir, gamma, mu, velocity_unit, nBox):
        '''
        Calculate the phase space for the a specific box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
            nBox (int): index of the box to load
        Returns:
            (arr): 2D histogram array of the overdensity-temperature phase space
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        assert self.Snap.head.DataHead.check_nbox(nBox)
        
        box_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure box_head has local_dims and offset
        assert box_head.head_set
        
        box = ChollaHydroBox(box_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(box_head.local_dims)
        
        # load density once as variable since we use it twice
        density = box.get_data(namebase, dataDir, 
                               self.Snap.head.DataHead.HydroHead.density_str, 
                               self.Snap.head.DataHead.old_format)
        
        return self.Analysis.create_phase(np.log10(box_calc.gas_temp(box.get_data(namebase, dataDir, 
                                                                                  self.Snap.head.DataHead.HydroHead.gasenergy_str, 
                                                                                  self.Snap.head.DataHead.old_format), 
                                                                     density, gamma, mu, velocity_unit)).ravel(),
                                          np.log10(box_calc.overdensity_median(density)).ravel()
                                         )
    
    
    def get_phasespace_loop(self, namebase, dataDir, gamma, mu, velocity_unit, nBoxes=None):
        '''
        Calculate the phase space for the entire sim box by looping over each
            individual box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D histogram array of the overdensity-temperature phase space
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        arr = self.Analysis.create_subarr((49,49))
        
        # grab requested hydro box heads
        hydro_box_heads = self.set_hboxheads(nBoxes)
        
        # grab density + gasenergy string
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
        
        for box_head in hydro_box_heads:
            # ensure box_head has local_dims and offset
            assert box_head.head_set
            
            # create a box + calculator object
            box = ChollaHydroBox(box_head, self.Snap.head.nSnap)
            box_calc = ChollaCalculator(box_head.local_dims)
            box_analysis = ChollaAnalysis(box_head.local_dims)
            
            # load density once as variable since we use it twice
            density = box.get_data(namebase, dataDir, 
                                   self.Snap.head.DataHead.HydroHead.density_str, 
                                   self.Snap.head.DataHead.old_format)
        
            phasespacebox, xedges, yedges = box_analysis.create_phase(np.log10(box_calc.gas_temp(box.get_data(namebase, dataDir, 
                                                                                                              self.Snap.head.DataHead.HydroHead.gasenergy_str, 
                                                                                                              self.Snap.head.DataHead.old_format), 
                                                                                                 density, gamma, mu, velocity_unit)).ravel(),
                                                                      np.log10(box_calc.overdensity_median(density)).ravel())
            
            arr += phasespacebox
        
        
        return arr, xedges, yedges
    
        
    def plot_cosmo_diagnostic_box(self, namebase, dataDir, gamma, mu, velocity_unit, nBox, show_ticks=True, fname=None):
        '''
        Plot the cosmological diagnostic for a specific box. Cosmological 
            diagnostic is defined in ChollaVizAnalysis
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            ...
        '''
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # make sure requested nBox is valid
        assert self.Snap.head.DataHead.check_nbox(nBox)
        
        box_head = self.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure box_head has local_dims and offset
        assert box_head.head_set
        
        # create a box + calculator + analysis + vizanalysis object
        box = ChollaHydroBox(box_head, self.Snap.head.nSnap)
        box_calc = ChollaCalculator(box_head.local_dims)
        box_analysis = ChollaAnalysis(box_head.local_dims)
        box_vizanalysis = ChollaVizAnalysis(box_head.local_dims, box_head.offset)
        
        # load density once as variable since we use it five times
        density = box.get_data(namebase, dataDir, 
                               self.Snap.head.DataHead.HydroHead.density_str, 
                               self.Snap.head.DataHead.old_format)
        
        hist, xedges, yedges = box_analysis.create_phase(np.log10(box_calc.gas_temp(box.get_data(namebase, dataDir, 
                                                                                                 self.Snap.head.DataHead.HydroHead.gasenergy_str, 
                                                                                                 self.Snap.head.DataHead.old_format), 
                                                                                    density, gamma, mu, velocity_unit)).ravel(), 
                                                         np.log10(box_calc.overdensity_median(density)).ravel())
        
        box_vizanalysis.cosmo_diagnostic(hist, xedges, yedges, 
                                         box_calc.densityk_projection(density, 0), 
                                         box_calc.densityk_projection(density, 1), 
                                         box_calc.densityk_projection(density, 2),
                                         show_ticks=show_ticks, fname=fname)
        
        
    def plot_cosmo_diagnostic_loop(self, namebase, dataDir, gamma, mu, velocity_unit, nBoxes = None, show_ticks=True, fname=None):
        '''
        Plot the cosmological diagnostic for a specific box. Cosmological 
            diagnostic is defined in ChollaVizAnalysis
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            ...
        '''
        
        assert self.Snap.head.DataHead.head_set
        assert self.Snap.head.DataHead.HydroHead.head_set
        
        # grab density + gasenergy string
        density_str = self.Snap.head.DataHead.HydroHead.density_str
        gasenergy_str = self.Snap.head.DataHead.HydroHead.gasenergy_str
        
        print("STARTING COSMO DIAGNOSTIC LOOP")
        time0 = time()
        
        dims = self.Snap.head.DataHead.HydroHead.dims
        n_x = self.Analysis.create_subarr((dims[1], dims[2]))
        n_y = self.Analysis.create_subarr((dims[0], dims[2]))
        n_z = self.Analysis.create_subarr((dims[0], dims[1]))
        
        # default phase space binning is (50,50)
        phasespace = self.Analysis.create_subarr((49,49))
        
        # grab requested hydro box heads
        hydro_box_heads = self.Snap.get_hydro_box_heads(nBoxes)
        
        time_boxes = np.zeros(self.Snap.head.DataHead.nBoxes)
        
        time_density = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_gasenergy = np.zeros(self.Snap.head.DataHead.nBoxes)
        
        time_overdensity = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_gastemp = np.zeros(self.Snap.head.DataHead.nBoxes)
        
        time_loggastemp = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_logoverdensity = np.zeros(self.Snap.head.DataHead.nBoxes)
        
        time_phase = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_nx = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_ny = np.zeros(self.Snap.head.DataHead.nBoxes)
        time_nz = np.zeros(self.Snap.head.DataHead.nBoxes)
        
        for i, hydro_box_head in enumerate(hydro_box_heads):
            print(f"Starting loop for box {i:.0f}")
            time_boxbox = time()
            
            # ensure hydro_box_head has local_dims and offset
            assert hydro_box_head.head_set
            box = ChollaHydroBox(hydro_box_head, self.Snap.head.nSnap)
            box_calc = ChollaCalculator(hydro_box_head.local_dims)
            box_analysis = ChollaAnalysis(hydro_box_head.local_dims)
            
            # load density once as variable since we use it five times
            time_box0 = time()
            density = box.get_data(namebase, dataDir, density_str, 
                                   self.Snap.head.DataHead.old_format)
            time_density[i] = time() - time_box0
            print(f"\tLoading density data took {time_density[i]:.2f} secs")
            
            # to test, look at how long loading gas energy takes
            time_box0 = time()
            gasenergy = box.get_data(namebase, dataDir, gasenergy_str, 
                                     self.Snap.head.DataHead.old_format)
            time_gasenergy[i] = time() - time_box0
            print(f"\tLoading gas energy data took {time_gasenergy[i]:.2f} secs")
            
            
            time_box0 = time()
            log_overdensity = box_calc.overdensity_median(density)
            time_overdensity[i] = time() - time_box0
            print(f"\tCalculating overdensity took {time_overdensity[i]:.2f} secs")
            
            time_box0 = time()
            log_gastemp = box_calc.gas_temp(gasenergy, density, gamma, mu, velocity_unit)
            time_gastemp[i] = time() - time_box0
            print(f"\tCalculating gas temp took {time_gastemp[i]:.2f} secs")
            
            
            time_box0 = time()
            log_overdensity = np.log10(log_overdensity).ravel()
            time_logoverdensity[i] = time() - time_box0
            print(f"\tTaking log of overdensity took {time_logoverdensity[i]:.2f} secs")
            
            time_box0 = time()
            log_gastemp = np.log10(log_gastemp).ravel()
            time_loggastemp[i] = time() - time_box0
            print(f"\tTaking log of gas temp took {time_loggastemp[i]:.2f} secs")
            
            
            # calculate phase space and projections
            time_box0 = time()
            phasespacebox, xedges, yedges = box_analysis.create_phase(log_gastemp, log_overdensity)
            time_phase[i] = time() - time_box0
            print(f"\tCalculating phase space took {time_phase[i]:.2f} secs")
            
            time_box0 = time()
            n_x_box = box_calc.densityk_projection(density, 0)
            time_nx[i] = time() - time_box0
            print(f"\tCalculating x-projection took {time_nx[i]:.2f} secs")
            
            time_box0 = time()
            n_y_box = box_calc.densityk_projection(density, 1)
            time_ny[i] = time() - time_box0
            print(f"\tCalculating y-projection took {time_ny[i]:.2f} secs")
            
            time_box0 = time()
            n_z_box = box_calc.densityk_projection(density, 2)
            time_nz[i] = time() - time_box0
            print(f"\tCalculating z-projection took {time_nz[i]:.2f} secs")
            
            # place box data onto overall dataset
            startX, startY, startZ = hydro_box_head.offset
            localX, localY, localZ = hydro_box_head.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            phasespace += phasespacebox
            n_x[startY:endY, startZ:endZ] += n_x_box
            n_y[startX:endX, startZ:endZ] += n_y_box
            n_z[startX:endX, startY:endY] += n_z_box
            
            time_boxes[i] = time() - time_boxbox
            print(f"\tEntire box analysis took {time_boxes[i]:.2f} secs")
        
        
        print(f"Each box density load took avg of {np.median(time_density):.2f} secs with std dev of {np.std(time_density):.2f} secs")
        print(f"Each box gas energy load took avg of {np.median(time_gasenergy):.2f} secs with std dev of {np.std(time_gasenergy):.2f} secs")
        
        print(f"Each box overdensity calculation took avg of {np.median(time_overdensity):.2f} secs with std dev of {np.std(time_overdensity):.2f} secs")
        print(f"Each box gas temp calculation took avg of {np.median(time_gastemp):.2f} secs with std dev of {np.std(time_gastemp):.2f} secs")
        
        print(f"Each box log overdensity took avg of {np.median(time_logoverdensity):.2f} secs with std dev of {np.std(time_logoverdensity):.2f} secs")
        print(f"Each box log gas temp took avg of {np.median(time_loggastemp):.2f} secs with std dev of {np.std(time_loggastemp):.2f} secs")
        
        print(f"Each box phase space calc took avg of {np.median(time_phase):.2f} secs with std dev of {np.std(time_phase):.2f} secs")
        print(f"Each box X-projection took avg of {np.median(time_nx):.2f} secs with std dev of {np.std(time_nx):.2f} secs")
        print(f"Each box Y-projection took avg of {np.median(time_ny):.2f} secs with std dev of {np.std(time_ny):.2f} secs")
        print(f"Each box Z-projection took avg of {np.median(time_nz):.2f} secs with std dev of {np.std(time_nz):.2f} secs")
        
        print(f"Each box took avg of {np.median(time_boxes):.2f} secs with std dev of {np.std(time_boxes):.2f} secs")
        
        # for entire box, offset will just be zeros
        chviz_analysis = ChollaVizAnalysis(self.Snap.head.DataHead.HydroHead.dims, (0.,0.,0.))
        
        print("calling cosmo diagnostic function")
        chviz_analysis.cosmo_diagnostic(phasespace, xedges, yedges, n_x, n_y, n_z,
                                        show_ticks=show_ticks, fname=fname)
        print(f"ENDING COSMO DIAGNOSTIC LOOP. Took total of {time() - time0:.2f} seconds")
