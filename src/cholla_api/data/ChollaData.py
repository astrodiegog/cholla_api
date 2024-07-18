import numpy as np
import h5py
from time import time

from cholla_api.data.ChollaHydroBox import ChollaHydroBoxHead
from cholla_api.data.ChollaParticleBox import ChollaParticleBoxHead

class ChollaParticleHead:
    '''
    Cholla Particle Data Head object
        Holds specific information regarding the particle data for a given 
            snapshot
        Initialized with:
        - nBoxes (int): total number of boxes used to run simulation
    
    Information about the particle box data children are saved in 
        ParticleBoxHeads, an ordered array of ChollaParticleBoxHead objects
    '''
    
    def __init__(self, nBoxes):
        self.ParticleBoxHeads = np.empty(nBoxes, dtype=object)
        self.head_set = False
        self.boxheads_set = False
        
    def set_head(self, nSnap, nBox, namebase, dataDir, old_format):
        '''
        Set the header attributes for this object.
        
        Args:
            nSnap (int): what snapshot number to use
            nBox (int): what box number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        fName = '{0}_particles.{1}.{2}'.format(nSnap, namebase, nBox)
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(nSnap) + '/' + fName
        fObj = h5py.File(fPath, 'r')
        self.dims = fObj.attrs['dims']
        fObj.close()
        
        self.n_parts_total = 0
        self.head_set = True
        
    def set_particleboxheads(self, nBoxes, nSnap, namebase, dataDir, old_format):
        '''
        Populate ParticleBoxHeads with ChollaParticleBoxHead objects.
        
        Args:
            nBoxes (int): total number of boxes used to run simulation
            nSnap (int): what snapshot number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        
        self.n_parts_total = 0
        time_boxes = np.zeros(nBoxes)

        for nBox in range(nBoxes):
            time_box1 = time()

            self.ParticleBoxHeads[nBox] = ChollaParticleBoxHead(nBox)
            self.ParticleBoxHeads[nBox].set_head(nSnap, nBox, namebase, dataDir, 
                                                 old_format)
            self.ParticleBoxHeads[nBox].offset_nparts = int(self.n_parts_total)
            if nSnap>0:
                # no particles in nsnap=0
                self.n_parts_total += self.ParticleBoxHeads[nBox].local_nparts

            time_boxes[nBox] = time() - time_box1
        
        print(f"Each ParticleBoxHead set took avg of {np.median(time_boxes):.2f} secs w std dev of {np.std(time_boxes):.2f} secs")

        self.boxheads_set = True
    
    def check_datakey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed
        
        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected particle keys
        '''
        
        # hard coded in https://github.com/cholla-hydro/cholla/wiki/Output
        parts_allkeys_set = {'density', 'particle_IDs', 'pos_x', 'pos_y', 
                             'pos_z', 'vel_x', 'vel_y', 'vel_z'}
        
        return data_key in parts_allkeys_set
        

class ChollaHydroHead:
    '''
    Cholla Hydro Data Head object
        Holds specific information regarding the hydro data for a given snapshot
        Initialized with:
        - nBoxes (int): total number of boxes used to run simulation
    
    Information about the hydro box data children are saved in 
        HydroBoxHeads, an ordered array of ChollaHydroBoxHead objects
    '''
    
    def __init__(self, nBoxes):
        self.HydroBoxHeads = np.empty(nBoxes, dtype=object)
        self.head_set = False
        self.boxheads_set = False
        self.energy_str = "Energy"
        self.density_str = "density"
        self.momx_str = "momentum_x"
        self.momy_str = "momentum_y"
        self.momz_str = "momentum_z"
        self.gasenergy_str = "GasEnergy"
        

    def set_head(self, nSnap, nBox, namebase, dataDir, old_format):
        '''
        Set the header attributes for this object.
        
        Args:
            nSnap (int): what snapshot number to use
            nBox (int): what box number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        #if nSnap == 0:
        #    # snapshot hydro data not well defined for snap=0
        #    return
        
        fName = '{0}.{1}.{2}'.format(nSnap, namebase, nBox)
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(nSnap) + '/' + fName
        fObj = h5py.File(fPath, 'r')
        self.dims = fObj.attrs['dims']
        fObj.close()
        
        self.head_set = True
        
    def set_hydroboxheads(self, nBoxes, nSnap, namebase, dataDir, old_format):
        '''
        Populate HydroBoxHeads with ChollaHydroBoxHead objects.
        
        Args:
            nBoxes (int): total number of boxes used to run simulation
            nSnap (int): what snapshot number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        
        time_boxes = np.zeros(nBoxes)
        for nBox in range(nBoxes):
            time_box1 = time()

            self.HydroBoxHeads[nBox] = ChollaHydroBoxHead(nBox)
            self.HydroBoxHeads[nBox].set_head(nSnap, namebase, dataDir, old_format)

            time_boxes[nBox] = time() - time_box1
        
        print(f"Each HydroBoxHead set took avg of {np.median(time_boxes):.2f} secs w std dev of {np.std(time_boxes):.2f} secs")
        self.boxheads_set = True        

    
    def check_datakey(self, data_key):
        '''
        Check if a requested data key is valid to be accessed
        
        Args:
            data_key (str): key string that will be used to access hdf5 dataset
        Return:
            (bool): whether data_key is a part of expected hydro keys
        '''
        
        # hard coded in https://github.com/cholla-hydro/cholla/wiki/Output
        hydro_allkeys_set = {"Energy", "density", "momentum_x", "momentum_y", 
                             "momentum_z", "GasEnergy"}
        
        return data_key in hydro_allkeys_set


        
        
class ChollaDataHead:
    '''
    Cholla Data Head object
        Holds specific information regarding the data within a snapshot
        Initialized with:
        - nBoxes (int): total number of boxes used to run simulation
        
    Information about the hydro and particle data children are saved in 
        HydroHead and ParticleHead
    '''
    
    def __init__(self, nBoxes, old_format=False):
        self.nBoxes = nBoxes
        self.old_format = old_format
        self.HydroHead = ChollaHydroHead(self.nBoxes)
        self.ParticleHead = ChollaParticleHead(self.nBoxes)
        self.head_set = False
        self.hydrohead_set = False
        self.particlehead_set = False
        
    def set_hydrohead(self, nSnap, namebase, dataDir):
        '''
        Set the HydroHead and populate its children HydroBoxHeads
        
        Args:
            nSnap (int): what snapshot number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            ...
        '''
        # use nBox=0 to set HydroHead
        self.HydroHead.set_head(nSnap, 0, namebase, dataDir, self.old_format)
        self.HydroHead.set_hydroboxheads(self.nBoxes, nSnap, namebase, dataDir, 
                                         self.old_format)
        self.hydrohead_set = True
        
    def set_particlehead(self, nSnap, namebase, dataDir):
        '''
        Set the ParticleHead and populate its children ParticleBoxHeads
        
        Args:
            nSnap (int): what snapshot number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            ...
        '''
        # use nBox=0 to set ParticleHead
        self.ParticleHead.set_head(nSnap, 0, namebase, dataDir, self.old_format)
        self.ParticleHead.set_particleboxheads(self.nBoxes, nSnap, namebase, 
                                               dataDir, self.old_format)
        self.particlehead_set = True
        
    def set_head(self, nSnap, namebase, dataDir, particles_flag):
        '''
        Set the HydroHead and (if flagged) ParticleHead.
        
        Args:
            nSnap (int): what snapshot number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            particles_flag (bool): whether particle data was saved
        Returns:
            ...
        '''
        
        self.set_hydrohead(nSnap, namebase, dataDir)
        if particles_flag:
            self.set_particlehead(nSnap, namebase, dataDir)
        self.head_set = True
            
    def check_nbox(self, nBox):
        '''
        Check if a requested nBox is valid to be indexed into.
        
        Args:
            nBox (int): index of the box to check
        Returns:
            (bool): whether 0 <= nBox < nBoxes
        '''
        
        if nBox >= self.nBoxes:
            print(f"{nBox:.0f} is an invalid box number. nBox must be less than {self.nBoxes:.0f} - the total number of boxes loaded onto DataHead.")
            return False
        elif nBox < 0:
            print(f"{nBox:.0f} is an invalid box number. nBox must be greater than or equal to 0.")
            return False
        return True
