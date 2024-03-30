import numpy as np
import h5py

from cholla_api.snap.ChollaSnap import ChollaSnapHead
from cholla_api.snap.ChollaSnap import ChollaSnap



class ChollaGlobalHead:
    '''
    Cholla Global Head object
        Holds global information regarding the entire simulation run.
        Initialized with:
        - dataPath (str): path to the data directory
        - namebase (str): middle string for data file names
        - nSnaps (int): total number of snapshots saved
    
    Information about the snapshot children are saved in SnapHeads, an ordered
        array of ChollaSnapHead objects
    '''
    
    def __init__(self, dataPath, namebase, nSnaps):
        self.dataDir = dataPath
        self.namebase = namebase
        self.nSnaps = nSnaps
        self.SnapHeads = np.empty(nSnaps, dtype=object)
        self.head_set = False
        self.snapheads_set = False
        
    def set_snapheads(self, nBoxes, particles_flag, cosmo_flag, old_format):
        '''
        Populate SnapHeads with ChollaSnapHead objects. Set each SnapHead's 
            object and its ChollaDataHead's object
        
        Args:
            nBoxes (int): total number of boxes used to run simulation
            particles_flag (bool): whether particle data was saved
            cosmo_flag (bool): whether cosmology type was used
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        for nSnap in range(self.nSnaps):
            self.SnapHeads[nSnap] = ChollaSnapHead(nSnap)
            if nSnap==0 and cosmo_flag:
                # file data not well defined for cosmo snap=0
                continue
            
            # use box=0 to set snaphead
            try:
                self.SnapHeads[nSnap].set_head(0, self.namebase, self.dataDir,
                                               particles_flag, cosmo_flag, old_format)
                try:
                    self.SnapHeads[nSnap].set_datahead(nBoxes, self.namebase, 
                                                   self.dataDir, particles_flag, old_format)
                except OSError:
                    print(f"Unable to open all {nBoxes:.0f} boxes for snapshot {nSnap:.0f}")
            except OSError:
                print(f"Unable to open snapshot {nSnap:.0f}. Setting SnapHead as None.")
                self.SnapHeads[nSnap] = None
            
            
        self.snapheads_set = True
    
    def set_head(self, nSnap=1, nBox=0, all_units=True, particles_flag=False, cosmo_flag=False, old_format=False):
        '''
        Set the header attributes for this object. Default use first new snapshot
            and first box.
        
        Args:
            nSnap (int): (optional) what snapshot number to use
            nBox (int): (optional) what box number to use
            all_units (bool): (optional) whether to include all units as attrs
            particles_flag (bool): (optional) whether particle data was saved
            cosmo_flag (bool): (optional) whether cosmology type was used
            old_format (bool): (optional) whether the old file structure was used
        Returns:
            ...
        '''
        fName = '{0}.{1}.{2}'.format(nSnap, self.namebase, nBox)
        if old_format:
            fPath = self.dataDir + '/' + fName
        else:
            fPath = self.dataDir + '/' + str(nSnap) + '/' + fName
        fObj = h5py.File(fPath, 'r')
                
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
            # Brant Robertson shared the velocity unit conversion
            self.velocity_unit = 1.0e5 # km to cm
        fObj.close()
        
        # grab particle info
        if particles_flag:
            fName = '{0}_particles.{1}.{2}'.format(nSnap, self.namebase, nBox)
            if old_format:
                fPath = self.dataDir + '/' + fName
            else:
                fPath = self.dataDir + '/' + str(nSnap) + '/' + fName
            fObj = h5py.File(fPath, 'r')
            self.paticle_mass_unit = fObj.attrs['particle_mass']
            fObj.close()
        
        self.head_set = True


class ChollaRun:
    '''
    Cholla Simulation Run object
        Holds information related to an entire simulation run. Initialized with
        - basePath (str): the directory the simulation is held in
        - nSnaps (int): total number of Snapshots saved
        - nBoxes (int): total number of boxes used to run simulation
        - namebase (str): (optional) middle string for the hdf5 data files
        - particles_flag (bool): (optional) whether particle data was saved
        - data_subdir (str): (optional) where data is placed within basePath
        - parts_flag (bool): (optional) whether particle data was saved
        - DE_flag (bool): (optional) whether dual energy formalism was used
        - cosmo_flag (bool): (optional) whether cosmology type was used
        - old_format (bool):(optional) whether the old file structure was used
        - nSnap_global (int): (optional) which snapshot to set global attributes
    
    '''

    def __init__(self, basePath, nSnaps, nBoxes, namebase='h5', data_subdir = '/data', parts_flag=False, DE_flag=False, cosmo_flag=False, old_format=False, nSnap_global=1):
        dataPath = basePath + data_subdir
        
        self.DE_flag = DE_flag
        self.parts_flag = parts_flag
        self.cosmo_flag = cosmo_flag
        self.old_format = old_format
        
        self.GlobalHead = ChollaGlobalHead(dataPath, namebase, nSnaps)
        self.GlobalHead.set_head(nSnap=nSnap_global, particles_flag=parts_flag, 
                                 cosmo_flag=cosmo_flag, old_format=old_format)
        
        self.GlobalHead.set_snapheads(nBoxes, parts_flag, cosmo_flag, old_format)
        
    def createSnap(self, nSnap):
        '''
        Lightweight method to create a ChollaSnap object with all keys from a 
            specific dataset
        
        Args:
            nSnap (int): the snapshot number to load
        Returns:
            ch_snap (ChollaSnap): snapshot object
        '''
        return ChollaSnap(self.GlobalHead.SnapHeads[nSnap])
        

