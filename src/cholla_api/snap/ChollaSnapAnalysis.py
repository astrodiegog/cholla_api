import numpy as np

from cholla_api.data.ChollaHydroBox import ChollaHydroBox

from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.snap.ChollaSnapCalc import ChollaSnapCalc
from cholla_api.analysis.ChollaAnalysis import ChollaAnalysis

class ChollaSnapAnalysis:
    '''
    Cholla Snapshot Analysis object
        This object will handle the arguments that are passed onto 
            ChollaAnalysis for some given SnapHead
        
        Similar to SnapCalc, where we will have global and nBox specific methods
    '''
    
    def __init__(self, SnapHead):
        self.SnapCalc = ChollaSnapCalc(SnapHead)
        
        self.Analysis = ChollaAnalysis(self.SnapCalc.Calc.dims)
    
    
    def set_hboxheads(self, nBoxes):
        '''
        Create the hydro box heads to loop over
        '''
        
        if nBoxes is None:
            return self.SnapCalc.Snap.head.DataHead.HydroHead.HydroBoxHeads
        else:
            hydro_box_heads = np.empty(len(nBoxes), dtype=object)
            for nBox, nBox_toload in enumerate(nBoxes):
                # make sure requested nBox is valid
                assert self.SnapCalc.Snap.head.DataHead.check_nbox(nBox_toload)
                
                hydro_box_heads[nBox] = self.SnapCalc.Snap.head.DataHead.HydroHead.HydroBoxHeads[nBox_toload]
            
            return hydro_box_heads
    
    
    def loop_gastemp_box(self, namebase, dataDir, gamma, mu, velocity_unit, nBoxes=None):
        '''
        Loop over nBoxes, apply gastemp_box, and place onto larger array using 
            the box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            gamma (float): ratio of specific heats
            mu (float): mean molecular weight in amu
            velocity_unit (float): unit conversion from velocity code units to cgs
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        arr = self.Analysis.create_arr()
        
        hydro_box_heads = self.set_hboxheads(nBoxes)
        
        for HydroBoxHead in hydro_box_heads:
            # ensure HydroBoxHead has local_dims and offset
            assert HydroBoxHead.head_set
            
            print(f"Starting with {HydroBoxHead.nBox:.0f}")
            data = self.SnapCalc.get_gastemp_box(namebase, dataDir, gamma, mu, 
                                                 velocity_unit, HydroBoxHead.nBox)
            
            ch_hydrobox = ChollaHydroBox(HydroBoxHead, self.SnapCalc.Snap.head.nSnap)
            
            ch_hydrobox.place_anydata(data, arr)
            print(f"Done with {HydroBoxHead.nBox:.0f}")
        
        return arr
    
    
    def loop_xproj_box(self, namebase, dataDir, nBoxes=None):
        '''
        Loop over nBoxes, apply x projection, and place onto larger array using 
            the box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        new_dims = (self.SnapCalc.Snap.head.DataHead.HydroHead.dims[1], 
                    self.SnapCalc.Snap.head.DataHead.HydroHead.dims[2])
        arr = self.Analysis.create_subarr(new_dims)
        
        hydro_box_heads = self.set_hboxheads(nBoxes)
        
        for HydroBoxHead in hydro_box_heads:
            # ensure HydroBoxHead has local_dims and offset
            assert HydroBoxHead.head_set
            
            data = self.SnapCalc.get_xprojection_box(namebase, dataDir, 
                                                     HydroBoxHead.nBox)
            
            startX, startY, startZ = HydroBoxHead.offset
            localX, localY, localZ = HydroBoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            arr[startY:endY, startZ:endZ] += data
        
        return arr
    
    def loop_yproj_box(self, namebase, dataDir, nBoxes=None):
        '''
        Loop over nBoxes, apply y projection, and place onto larger array using 
            the box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        new_dims = (self.SnapCalc.Snap.head.DataHead.HydroHead.dims[0], 
                    self.SnapCalc.Snap.head.DataHead.HydroHead.dims[2])
        arr = self.Analysis.create_subarr(new_dims)
        
        hydro_box_heads = self.set_hboxheads(nBoxes)
        
        for HydroBoxHead in hydro_box_heads:
            # ensure HydroBoxHead has local_dims and offset
            assert HydroBoxHead.head_set
            
            data = self.SnapCalc.get_yprojection_box(namebase, dataDir, 
                                                     HydroBoxHead.nBox)
            
            startX, startY, startZ = HydroBoxHead.offset
            localX, localY, localZ = HydroBoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            arr[startX:endX, startZ:endZ] += data
        
        return arr
    
    
    def loop_zproj_box(self, namebase, dataDir, nBoxes=None):
        '''
        Loop over nBoxes, apply z projection, and place onto larger array using 
            the box's offset
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            ...
        '''
        
        new_dims = (self.SnapCalc.Snap.head.DataHead.HydroHead.dims[0], 
                    self.SnapCalc.Snap.head.DataHead.HydroHead.dims[1])
        arr = self.Analysis.create_subarr(new_dims)
        
        hydro_box_heads = self.set_hboxheads(nBoxes)
        
        for HydroBoxHead in hydro_box_heads:
            # ensure HydroBoxHead has local_dims and offset
            assert HydroBoxHead.head_set
            
            data = self.SnapCalc.get_zprojection_box(namebase, dataDir, 
                                                     HydroBoxHead.nBox)
            
            startX, startY, startZ = HydroBoxHead.offset
            localX, localY, localZ = HydroBoxHead.local_dims
            endX, endY, endZ = startX + localX, startY + localY, startZ + localZ
            
            arr[startX:endX, startY:endY] += data
        
        return arr
    
    
    def get_phasespace(self, namebase, dataDir, gamma, mu, velocity_unit, nBoxes=None):
        '''
        Calculate the phase space for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            (arr): 2D histogram array of the overdensity-temperature phase space
        '''
        
        return self.Analysis.create_phase(np.log10(self.SnapCalc.get_gastemp(namebase, dataDir, 
                                                                gamma, mu, velocity_unit, 
                                                                nBoxes).flatten()),
                                          np.log10(self.SnapCalc.get_overdensity(namebase, 
                                                                dataDir, nBoxes).flatten()))
    
    
    def get_phasespace_box(self, namebase, dataDir, gamma, mu, velocity_unit, nBox):
        '''
        Calculate the phase space for the entire sim box
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            nBox (int): index of the box to load
        Returns:
            (arr): 2D histogram array of the overdensity-temperature phase space
        '''
        
        return self.Analysis.create_phase(np.log10(self.SnapCalc.get_gastemp_box(namebase, dataDir, 
                                                                gamma, mu, velocity_unit, 
                                                                nBox).flatten()),
                                          np.log10(self.SnapCalc.get_overdensity_box(namebase,
                                                                dataDir, nBox).flatten()))
