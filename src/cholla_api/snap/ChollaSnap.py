import numpy as np
import h5py

from cholla_api.data.ChollaData import ChollaDataHead
from cholla_api.data.ChollaHydroBox import ChollaHydroBox
from cholla_api.data.ChollaParticleBox import ChollaParticleBox

class ChollaSnapHead:
    '''
    Cholla Snapshot Head object
        Holds snapshot specific information
        Initialized with:
        - nSnap (int): number of the snapshot within run
    '''
    
    def __init__(self, nSnap):
        self.nSnap = nSnap
        self.head_set = False
        self.datahead_set = False
    
    def set_head(self, nBox, namebase, dataDir, particles_flag, cosmo_flag, old_format):
        '''
        Set the header attributes for this object
        
        Args:
            nBox (int): what box number to use
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            particles_flag (bool): whether particle data was saved
            cosmo_flag (bool): whether cosmology type was used
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        fName = '{0}.{1}.{2}'.format(self.nSnap, namebase, nBox)
        if old_format:
            fPath = dataDir + '/' + fName
        else:
            fPath = dataDir + '/' + str(self.nSnap) + '/' + fName
        fObj = h5py.File(fPath, 'r')
        
        # grab t + dt
        self.t = float(fObj.attrs['t'])
        self.dt = float(fObj.attrs['dt'])
        
        # grab cosmology info
        # cosmo info not well defined for snap=0
        if cosmo_flag and self.nSnap>0:
            self.a = float(fObj.attrs['Current_a'])
            self.z = float(fObj.attrs['Current_z'])
        
        fObj.close()
        
        # grab particles info
        # particle info not well defined for snap=0
        if particles_flag and self.nSnap>0:
            fName = '{0}_particles.{1}.{2}'.format(self.nSnap, namebase, nBox)
            if old_format:
                fPath = dataDir + '/' + fName
            else:
                fPath = dataDir + '/' + str(self.nSnap) + '/' + fName
            fObj = h5py.File(fPath, 'r')
            self.t_particles = fObj.attrs['t_particles']
            self.dt_particles = fObj.attrs['dt_particles']
            fObj.close()
        
        self.head_set = True
        
    
    def set_datahead(self, nBoxes, namebase, dataDir, particles_flag, old_format):
        '''
        Create and populate DataHead object. Set each DataHydroHeads's and
            DataParticleHeads's arrays of objects
        
        Args:
            nBoxes (int): total number of boxes used to run simulation
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            particles_flag (bool): whether particle data was saved
            old_format (bool): whether the old file structure was used
        Returns:
            ...
        '''
        self.DataHead = ChollaDataHead(nBoxes, old_format)
        self.DataHead.set_head(self.nSnap, namebase, dataDir, particles_flag)
        self.datahead_set = True



class ChollaSnap:
    '''
    Cholla Snapshot object
        Holds snapshot specific information and methods to access data for that
            snapshot
        
        My two main get_data() functions loop over the BoxHeads and place all
            of the requested data onto one array.
    
    TODO:
        create method to load data within a requested domain. would I have to 
            create a separate method in BoxHeads for global domain? this would
            save having to save the local dims and offset redundant.
        LOOK AT BRUNO'S PULL REQUEST. he essentially does something very similar
            to this!
    '''
    
    def __init__(self, SnapHead):
        self.head = SnapHead
        # assert that the SnapHead is already set
        assert self.head.head_set
        
    def get_hydroboxdata(self, namebase, dataDir, data_key, nBox):
        '''
        Grab data for a specific Box
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            nBox (int): index of the box to load
        Returns:
            (arr): array that will hold data information
        '''
        assert self.head.DataHead.head_set
        assert self.head.DataHead.HydroHead.head_set
        
        # make sure requested nBox is valid
        assert self.head.DataHead.check_nbox(nBox)
        
        # make sure requested data_key is valid
        assert self.head.DataHead.HydroHead.check_datakey(data_key)
        
        HydroBoxHead = self.head.DataHead.HydroHead.HydroBoxHeads[nBox]
        
        # ensure HydroBoxHead has local_dims and offset
        assert HydroBoxHead.head_set
        ch_hydrobox = ChollaHydroBox(HydroBoxHead, self.head.nSnap)
        
        return ch_hydrobox.get_data(namebase, dataDir, data_key, self.head.DataHead.old_format)
    
    
    def get_particleboxdata(self, namebase, dataDir, data_key, nBox):
        '''
        Grab data for a specific Box
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            nBox (int): index of the box to load
        Returns:
            arr (arr): array that will hold data information
        '''
        assert self.head.DataHead.head_set
        assert self.head.DataHead.ParticleHead.head_set
        
        # make sure requested nBox is valid
        assert self.head.DataHead.check_nbox(nBox)
        
        # make sure requested data_key is valid
        assert self.head.DataHead.ParticleHead.check_datakey(data_key)
        
        ParticleBoxHead = self.head.DataHead.ParticleHead.ParticleBoxHeads[nBox]
        
        # ensure ParticleBoxHead has local_dims and offset
        assert ParticleBoxHead.head_set
        ch_particlebox = ChollaParticleBox(ParticleBoxHead, self.head.nSnap)
        
        return ch_particlebox.get_data(namebase, dataDir, data_key, self.head.DataHead.old_format)
    
    
    def get_hydrodata(self, namebase, dataDir, data_key, nBoxes=None):
        '''
        Grab and concatenate data onto one large array for the entire simulation
            box.
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            arr (arr): array that will hold data information
        '''
        assert self.head.DataHead.head_set
        assert self.head.DataHead.HydroHead.head_set
        
        # make sure requested data_key is valid
        assert self.head.DataHead.HydroHead.check_datakey(data_key)
        
        # initialize array
        arr = np.zeros(self.head.DataHead.HydroHead.dims)
        
        # create the hydro box heads to loop over
        if nBoxes is not None:
            hydro_box_heads = np.empty(len(nBoxes), dtype=object)
            for nBox, nBox_toload in enumerate(nBoxes):
                # make sure requested nBox is valid
                assert self.head.DataHead.check_nbox(nBox_toload)
                
                hydro_box_heads[nBox] = self.head.DataHead.HydroHead.HydroBoxHeads[nBox_toload]
        else:
            hydro_box_heads = self.head.DataHead.HydroHead.HydroBoxHeads
        
        for HydroBoxHead in hydro_box_heads:
            # ensure HydroBoxHead has local_dims and offset
            assert HydroBoxHead.head_set
            
            ch_hydrobox = ChollaHydroBox(HydroBoxHead, self.head.nSnap)
            ch_hydrobox.place_data(namebase, dataDir, data_key, arr,
                                   self.head.DataHead.old_format)
        
        return arr

    def get_particledata(self, namebase, dataDir, data_key, nBoxes=None):
        '''
        Grab and concatenate data onto one large array for the entire simulation
            box.
            
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
            data_key (str): key to access data from hdf5 file
            nBoxes (list): (optional) list of boxes to load, default is all
        Returns:
            arr (arr): array that will hold data information
        '''
        assert self.head.DataHead.head_set
        assert self.head.DataHead.ParticleHead.head_set
        
        # make sure requested data_key is valid
        assert self.head.DataHead.ParticleHead.check_datakey(data_key)
        
        densityCIC_key = "density"
        # initialize array
        if data_key == densityCIC_key:
            arr = np.zeros(self.head.DataHead.ParticleHead.dims)
        else:
            arr = np.zeros(self.head.DataHead.ParticleHead.n_parts_total)
        
        # create the particle box heads to loop over
        if nBoxes is not None:
            particle_box_heads = np.empty(len(nBoxes), dtype=object)
            for nBox, nBox_toload in enumerate(nBoxes):
                # make sure requested nBox is valid
                assert self.head.DataHead.check_nbox(nBox_toload)
                
                particle_box_heads[nBox] = self.head.DataHead.ParticleHead.ParticleBoxHeads[nBox_toload]
        else:
            particle_box_heads = self.head.DataHead.ParticleHead.ParticleBoxHeads
        
        for ParticleBoxHead in particle_box_heads:
            # ensure ParticleBoxHead has local_dims and offset
            assert ParticleBoxHead.head_set
            
            ch_particlebox = ChollaParticleBox(ParticleBoxHead, self.head.nSnap)
            ch_particlebox.place_data(namebase, dataDir, data_key, arr,
                                      self.head.DataHead.old_format)
        
        return arr
    
    def get_densityCIC(self, namebase, dataDir):
        '''
        Get particle density cloud-in-cell data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold density information
        '''
        
        return self.get_particledata(namebase, dataDir, "density")
    
    
    def get_posx(self, namebase, dataDir):
        '''
        Get particle x position data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold x position information
        '''
        
        return self.get_particledata(namebase, dataDir, "pos_x")
    
    
    def get_posy(self, namebase, dataDir):
        '''
        Get particle y position data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold x position information
        '''
        
        return self.get_particledata(namebase, dataDir, "pos_y")
    
    def get_posz(self, namebase, dataDir):
        '''
        Get particle z position data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold z position information
        '''
        
        return self.get_particledata(namebase, dataDir, "pos_z")
    
    
    def get_velx(self, namebase, dataDir):
        '''
        Get particle x velocity data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold x velocity information
        '''
        
        return self.get_particledata(namebase, dataDir, "vel_x")
    
    
    def get_vely(self, namebase, dataDir):
        '''
        Get particle y velocity data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold y velocity information
        '''
        
        return self.get_particledata(namebase, dataDir, "vel_y")
    
    def get_velz(self, namebase, dataDir):
        '''
        Get particle z velocity data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold z velocity information
        '''
        
        return self.get_particledata(namebase, dataDir, "vel_z")
    
    def get_partid(self, namebase, dataDir):
        '''
        Get particle id data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold particle id information
        '''
        
        return self.get_particledata(namebase, dataDir, "particle_IDs")
    
    
    def get_density(self, namebase, dataDir):
        '''
        Get hydro density data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold density information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "density")
    
    
    def get_momentumx(self, namebase, dataDir):
        '''
        Get hydro momentum in x-direction data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold momentum in x-direction information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "momentum_x")
    
    def get_momentumy(self, namebase, dataDir):
        '''
        Get hydro momentum in y-direction data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold momentum in y-direction information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "momentum_y")
    
    def get_momentumz(self, namebase, dataDir):
        '''
        Get hydro momentum in z-direction data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold momentum in z-direction information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "momentum_z")
    
    def get_energy(self, namebase, dataDir):
        '''
        Get hydro energy data
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold energy information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "Energy")
    
    
    def get_gasenergy(self, namebase, dataDir):
        '''
        Get hydro gas energy data. Only works if Dual Energy formalism was used
        
        Args:
            namebase (str): middle string for data file names
            dataDir (str): path to the data directory
        Returns:
            (arr): array that will hold gas energy information
        '''
        
        return self.get_hydrodata(namebase, dataDir, "GasEnergy")
    