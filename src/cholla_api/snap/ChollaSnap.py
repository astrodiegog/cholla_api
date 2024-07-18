import numpy as np
import h5py

from cholla_api.data.ChollaBox import ChollaBox

class ChollaSnapHead:
    '''
    Cholla Snapshot Head object
        Holds snapshot specific information
        Initialized with:
        - nSnap (int): number of the snapshot within run
    '''
    
    def __init__(self, nSnap):
        self.nSnap = nSnap
   
    def set_timeinfo(self, ChollaBox):
        '''
        Set time attributes for this object

        Args:
            ChollaBox (ChollaBox): what box to use to set attributes
        Returns:
            ...
        '''

        fPath = ChollaBox.get_hydrofPath()
        fObj = h5py.File(fPath, 'r')

        self.t = fObj.attrs['t'].item()
        self.dt = fObj.attrs['dt'].item()

        fObj.close()

    def set_cosmoinfo(self, ChollaBox):
        '''
        Set cosmology time attributes for this object

        Args:
            ChollaBox (ChollaBox): what box to use to set attributes
        Returns:
            ...
        '''

        fPath = ChollaBox.get_hydrofPath()
        fObj = h5py.File(fPath, 'r')

        self.a = fObj.attrs['Current_a'].item()
        self.z = fObj.attrs['Current_z'].item()

        fObj.close()

    def set_particleinfo(self, ChollaBox):
        '''
        Set particle time attributes for this object

        Args:
            ChollaBox (ChollaBox): what box to use to set attributes
        Returns:
            ...
        '''

        fPath = ChollaBox.get_particlefPath()
        fObj = h5py.File(fPath, 'r')

        self.t_particles = fObj.attrs['t_particles'].item()
        self.dt_particles = fObj.attrs['dt_particles'].item()

        fObj.close()



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
    '''
    
    def __init__(self, RunPath, ChollaSnapHead):
        self.SnapHead = ChollaSnapHead
        self.SnapPath = RunPath + '/' + str(self.SnapHead.nSnap)

    def get_hydroboxdata(self, ChollaBoxHead, key, dtype=np.float32):
        '''
        Grab hydro data for a specific Box

        Args:
            ChollaBoxHead (ChollaBoxHead): head for box to use
            key (str): key to access data from hydro hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns
            (arr): array holding data
        '''

        box = ChollaBox(self.SnapPath, ChollaBoxHead)

        return box.get_hydrodata(key, dtype)

    def get_particleboxdata(self, ChollaBoxHead, key, dtype=np.float32):
        '''
        Grab particle data for a specific Box

        Args:
            ChollaBoxHead (ChollaBoxHead): head for box to use
            key (str): key to access data from particle hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns
            (arr): array holding data
        '''

        box = ChollaBox(self.SnapPath, ChollaBoxHead)

        return box.get_particledata(key, dtype)

    def get_hydrodata(self, ChollaGrid, key, dtype=np.float32):
        '''
        Grab hydro data for all boxes and concatenate data

        Args:
            ChollaGrid (ChollaGrid): grid holding boxheads
            key (str): key to access data from hydro hdf5 file
        Returns:
            arr (arr): array holding global data
        '''

        arr = np.zeros((ChollaGrid.nx_global, ChollaGrid.ny_global, ChollaGrid.nz_global), dtype=dtype)

        for boxhead in ChollaGrid.get_BoxHeads():
            box = ChollaBox(self.SnapPath, boxhead)
            boxdata = box.get_hydrodata(key, dtype)
            box.place_data(boxdata, arr)

        return arr

    def get_particledata(self, ChollaGrid, key, dtype=np.float32):
        '''
        Grab particle data for all boxes and concatenate data

        Args:
            ChollaGrid (ChollaGrid): grid holding boxheads
            key (str): key to access data from particle hdf5 file
        Returns:
            arr (arr): array holding global data
        '''

        if (key == "density"):
            arr = np.zeros((ChollaGrid.nx_global, ChollaGrid.ny_global, ChollaGrid.nz_global), dtype=dtype)

            for boxhead in ChollaGrid.get_BoxHeads():
                box = ChollaBox(self.SnapPath, boxhead)
                boxdata = box.get_particledata(key, dtype)
                box.place_data(boxdata, arr)
        else:
            nparts_tot = 0
            nparts_local = np.zeros(ChollaGrid.nprocs, dtype=int)
            nparts_offsets = np.zeros(ChollaGrid.nprocs, dtype=int)

            for boxhead in ChollaGrid.get_BoxHeads():
                box = ChollaBox(self.SnapPath, boxhead)
                nparts_offsets[box.BoxHead.nBox] = nparts_tot
                nparts_local[box.BoxHead.nBox] = box.get_nparts()
                nparts_tot += nparts_local[box.BoxHead.nBox]

            arr = np.zeros(int(nparts_tot), dtype=dtype)

            for boxhead in ChollaGrid.get_BoxHeads():
                box = ChollaBox(self.SnapPath, boxhead)
                nStart = nparts_offsets[box.BoxHead.nBox]
                nLocal = nparts_local[box.BoxHead.nBox]
                nEnd = nStart + nLocal
                arr[nStart:nEnd] = box.get_particle_data(key, dtype)

        return arr
        
    
    
