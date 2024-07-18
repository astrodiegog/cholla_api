import numpy as np
from time import time
import h5py

from cholla_api.snap.ChollaSnap import ChollaSnapHead
from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.data.ChollaBox import ChollaBox

a = "In the future, I can just hold a MacroFlags"


class ChollaGlobal:
    '''
    Cholla Simulation Run object
        Holds information related to an entire simulation run. Initialized with
        - basePath (str): the directory the simulation is held in
        - chGrid (ChollaGrid): ChollaGrid object, holding domain information
        - gravity_flag (bool): whether gravity data was saved
        - particles_flag (bool): whether particle data was saved
        - cosmo_flag (bool): (optional) whether cosmology type was used
        - data_subdir (str): (optional) where data is placed within basePath
    '''

    def __init__(self, basePath,  ChollaGrid, gravity_flag, particles_flag, cosmo_flag, data_subdir = '/data'):

        self.dataPath = basePath + data_subdir
        
        self.grid = ChollaGrid
        self.gravity_flag = gravity_flag
        self.particles_flag = particles_flag
        self.cosmo_flag = cosmo_flag

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

        return ChollaBox(snap.SnapPath, boxhead)

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
        snap.SnapHead.set_time_info(box)

        if self.cosmo_flag:
            snap.SnapHead.set_cosmoinfo(box)
        if self.particles_flag:
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
        self.length_unit = fObj.attrs['length_unit']
        self.mass_unit = fObj.attrs['mass_unit']
        self.time_unit = fObj.attrs['time_unit']
        # following values can be reconstructed from previous base units
        # but can include if wanted :P
        if all_units:
            self.density_unit = fObj.attrs['density_unit']
            self.energy_unit = fObj.attrs['energy_unit']
            self.velocity_unit = fObj.attrs['velocity_unit']

        # grab dimensions and domain
        self.dims = fObj.attrs['dims']
        self.domain = fObj.attrs['domain']

        # grab dx + gamma
        self.dx = fObj.attrs['dx']
        self.gamma = fObj.attrs['gamma']

        # grab cosmology params
        if cosmo_flag:
            self.H0 = fObj.attrs['H0']
            self.Omega_L = fObj.attrs['Omega_L']
            self.Omega_M = fObj.attrs['Omega_M']
        
        fObj.close()

        # grab particle info
        if particles_flag:
            fObj = h5py.File(box.get_particlefPath(), 'r')
            self.paticle_mass_unit = fObj.attrs['particle_mass']
            fObj.close()


