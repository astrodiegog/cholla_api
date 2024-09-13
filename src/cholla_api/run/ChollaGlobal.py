import numpy as np
import h5py

from cholla_api.run.ChollaMacroFlags import ChollaMacroFlags
from cholla_api.snap.ChollaSnap import ChollaSnapHead
from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.data.ChollaBox import ChollaBox
from cholla_api.OTFanalysis.ChollaOnTheFlyAnalysis import ChollaOnTheFlyAnalysis
from cholla_api.OTFanalysis.ChollaOnTheFlySkewers import ChollaOnTheFlySkewers


class ChollaGlobal:
    '''
    Cholla Simulation Run object
        Holds information related to an entire simulation run. Initialized with
        - basePath (str): the directory the simulation is held in
        - chGrid (ChollaGrid): ChollaGrid object, holding domain information
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro compiling information
        - data_subdir (str): (optional) where data is placed within basePath
    '''

    def __init__(self, basePath,  chGrid, chMacroFlags, data_subdir = '/data', analysis_subdir = '/analysis', skewers_subdir = '/skewers'):
        self.dataPath = basePath + data_subdir
        self.grid = chGrid
        self.macroFlags = chMacroFlags
        if self.macroFlags.OTF_Analysis:
            self.analysisPath = basePath + analysis_subdir
            self.skewersPath = basePath + skewers_subdir


    def test_domaindecomp(self, nSnap=1):
        '''
        Make sure the ChollaGrid domain decomposition agrees with hydro files

        Args:
            nSnap (int): (optional) snapshot number to compare against
        Returns:
            (bool): whether a local domain & offset combination doesn't match up
        '''
        
        snaphead = ChollaSnapHead(nSnap)
        snap = ChollaSnap(self.dataPath, snaphead)

        for boxhead in self.grid.get_BoxHeads():
            box = ChollaBox(snap.SnapPath, boxhead)
            fObj = h5py.File(box.get_hydrofPath(), 'r')
            file_localdim = fObj.attrs['dims_local']
            file_offset = fObj.attrs['offset']
            fObj.close()

            if (not np.array_equal(file_localdim, boxhead.local_dims)):
                print(f"Failed local dim at boxhead {boxhead.nBox:.0f}")
                print(f"\t File has {file_localdim}")
                print(f"\t Box has {boxhead.local_dims}")
                return False
            if (not np.array_equal(file_offset, boxhead.offset)):
                print(f"Failed at offset boxhead {boxhead.nBox:.0f}")
                print(f"\t File has {file_offset}")
                print(f"\t Box has {boxhead.offset}")
                return False

        return True

    def get_cosmoHeader(self, nAnalysisOutput=0):
        '''
        Given an analysis output number, create the ChollaOnTheFlyAnalysis
            object to return the Cosmology Header object

        Args:
            nAnalysisOutput (int): (optional) analysis output
        Returns:
            (ChollaCosmologyHead): cosmology header object
        '''
        OTFAnalysis = self.get_OTFAnalysis(nAnalysisOutput)

        return OTFAnalysis.get_cosmoHeader()


    def get_OTFAnalysis(self, nAnalysisOutput):
        '''
        Given an analysis output number, return the ChollaOnTheFlyAnalysis
            object associated with the analysis subdirectory in initialization

        Args:
            nAnalysisOutput (int): analysis output
        Returns:
            (ChollaOnTheFlyAnalysis): OTF Analysis object
        '''

        assert self.macroFlags.OTF_Analysis

        return ChollaOnTheFlyAnalysis(nAnalysisOutput, self.analysisPath, self.grid)


    def get_OTFSkewers(self, nSkewersOutput):
        '''
        Given a skewer output number, return the ChollaOnTheFlySkewers
            object associated with the skewer subdirectory in initialization

        Args:
            nSkewersOutput (int): skewer output
        Returns:
            (ChollaOnTheFlySkewers): OTF Skewers object
        '''

        assert self.macroFlags.OTF_Analysis

        return ChollaOnTheFlySkewers(nSkewersOutput, self.skewersPath, self.grid)


    def get_snapbox(self, nSnap, nBox):
        '''
        Given a snapshot number and box number, return the ChollaBox object

        Args:
            nSnap (int): snapshot number to load
            nBox (int): process / box number to load
        Returns:
            (ChollaBox): box object
        '''
        # ensure nBox is valid
        assert nBox < self.grid.nprocs

        snaphead = ChollaSnapHead(nSnap)
        snap = ChollaSnap(self.dataPath, snaphead)
        boxhead = self.grid.get_BoxHead(nBox)

        return ChollaBox(snap.SnapPath, boxhead, self.macroFlags)

    def get_snap(self, nSnap, box_for_attrs=0):
        '''
        Given snapshot number, return the ChollaSnap object

        Args:
            nSnap (int): snapshot number to load
            box_for_attrs (int): (optional) box number to use for snapshot attributes
        Returns:
            snap (ChollaSnap): snapshot object
        '''

        snaphead = ChollaSnapHead(nSnap)
        snap = ChollaSnap(self.dataPath, snaphead)
        box = self.get_snapbox(nSnap, box_for_attrs)
        snap.SnapHead.set_timeinfo(box)

        if self.macroFlags.Cosmology:
            snap.SnapHead.set_cosmoinfo(box)
        if self.macroFlags.Particles:
            snap.SnapHead.set_particleinfo(box)

        return snap

    def set_globalattrs(self, nSnap, nBox, all_units=False):
        '''
         Args:
            nSnap (int): snapshot number to use for attributes
            nBox (int): process / box number to use for attributes
            all_units (bool): (optional) whether to save all units as attributes
        Returns:
            ...
        '''

        box = self.get_snapbox(nSnap, nBox)

        fObj = h5py.File(box.get_hydrofPath(), 'r')

        # unit values used to convert from code units to cgs values
        self.length_unit = fObj.attrs['length_unit'].item()
        self.mass_unit = fObj.attrs['mass_unit'].item()
        self.time_unit = fObj.attrs['time_unit'].item()
        # following values can be reconstructed from previous base units
        # but can include if wanted :P
        if all_units:
            self.density_unit = fObj.attrs['density_unit'].item()
            self.energy_unit = fObj.attrs['energy_unit'].item()
            self.velocity_unit = fObj.attrs['velocity_unit'].item()

        # grab dimensions and domain
        self.dims = fObj.attrs['dims']
        self.domain = fObj.attrs['domain']

        # grab dx + gamma
        self.dx = fObj.attrs['dx']
        self.gamma = fObj.attrs['gamma'].item()

        # grab cosmology params
        if self.macroFlags.Cosmology:
            self.H0 = fObj.attrs['H0'].item()
            self.Omega_L = fObj.attrs['Omega_L'].item()
            self.Omega_M = fObj.attrs['Omega_M'].item()
        
        fObj.close()

        # grab particle info
        if self.macroFlags.Particles:
            fObj = h5py.File(box.get_particlefPath(), 'r')
            self.paticle_mass_unit = fObj.attrs['particle_mass'].item()
            fObj.close()


