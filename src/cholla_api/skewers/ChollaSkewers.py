import numpy as np
import h5py

from cholla_api.data.ChollaBox import ChollaBox

class ChollaSkewerLocalFaceHead:
    '''
    Cholla Skewer Local Face Head

    Holds information regarding a local skewer within face

        Initialized with:
        - localface_id (int): id of the local skewer on face
        - localface_joffset (int): offset along j-axis
        - localface_koffset (int): offset along k-axis
    '''
    def __init__(self, localface_id, localface_joffset, localface_koffset):
        self.localface_id = localface_id
        self.localface_joffset = localface_joffset
        self.localface_koffset = localface_koffset


class ChollaSkewerFaceHead:
    '''
    Cholla Skewer Face Head

    Holds information regarding a local skewer face

        Initialized with:
        - face_id (int): id of the face
        - face_joffset (int): offset along j-axis
        - face_koffset (int): offset along k-axis
    '''
    def __init__(self, face_id, face_joffset, face_koffset):
        self.face_id = face_id
        self.face_joffset = face_joffset
        self.face_koffset = face_koffset


class ChollaSkewerGlobalHead:
    '''
    Cholla Skewer Global Head

    Holds information regarding a global skewer

        Initialized with:
        - global_id (int): id of the global skewer
        - chFaceHead (ChollaSkewerFaceHead): ChollaSkewerFaceHead object,
            holds info on face within grid
        - chFaceLocalHead (ChollaSkewerLocalFaceHead): ChollaSkewerLocalFaceHead
            object, holds info on skewer within face
        - n_los (int): number of cells along line-of-sight
        - nlos_proc (int): number of processes along line-of-sight
    '''
    def __init__(self, global_id, ChollaSkewerFaceHead, ChollaSkewerLocalFaceHead, n_los, nlos_proc):
        self.global_id = global_id
        self.skewFaceHead = ChollaSkewerFaceHead
        self.skewLocalFaceHead = ChollaSkewerLocalFaceHead
        self.n_los = int(n_los)
        self.nlos_proc = int(nlos_proc)

    def get_globalj(self):
        '''
        Grab global j offset

        Args:
            ...
        Returns:
            (int): global j offset
        '''

        return int(self.skewFaceHead.face_joffset + self.skewLocalFaceHead.localface_joffset)

    def get_globalk(self):
        '''
        Grab global k offset

        Args:
            ...
        Returns:
            (int): global k offset
        '''

        return int(self.skewFaceHead.face_koffset + self.skewLocalFaceHead.localface_koffset)



class ChollaSkewerxGlobal:
    '''
    Cholla Skewer x Global object
        Holds information regarding skewer along x direction including methods
            to grab data along skewer

        Initialized with:
        - ChollaSkewerGlobalHead (ChollaSkewerGlobalHead): global head object
        - SnapPath (str): path to snapshot directory
        - ChollaGrid (ChollaGrid): grid object
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, ChollaGrid):
        self.skewGlobalHead = ChollaSkewerGlobalHead
        self.SnapPath = SnapPath

        self.nx_local = int(self.skewGlobalHead.n_los / self.skewGlobalHead.nlos_proc)

        # save all of the ChollaBoxHead objects that this skewer passes through
        boxheads = self.skewGlobalHead.nlos_proc * [None]

        # grab j, k coordinates
        j = self.skewGlobalHead.get_globalj()
        k = self.skewGlobalHead.get_globalk()

        for curr_nlosproc in range(self.skewGlobalHead.nlos_proc):
            # starting x coordinate for skewer within current LOS process
            xstart = int(self.nx_local * curr_nlosproc)
            boxhead = ChollaGrid.get_BoxHead_ijk(xstart, j, k)
            boxheads[curr_nlosproc] = boxhead
        
        self.boxheads = tuple(boxheads)


    def get_hydrodata(self, key, dtype=np.float32):
        '''
        Return a specific hydro dataset along skewer

        Args:
            key (str): key to access data from hydro hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        arr = np.zeros(self.skewGlobalHead.n_los, dtype=dtype)
        j_local = self.skewGlobalHead.skewLocalFaceHead.localface_joffset
        k_local = self.skewGlobalHead.skewLocalFaceHead.localface_koffset

        for boxhead in self.boxheads:
            box = ChollaBox(self.SnapPath, boxhead)

            # check key with boxhead having 0 x-offset
            if (boxhead.offset[0] == 0):
                assert box.check_hydrokey(key)

            fObj = h5py.File(box.get_hydrofPath(), 'r')
            arr[boxhead.offset[0] : boxhead.offset[0] + boxhead.local_dims[0]] = fObj.get(key)[:, j_local, k_local]
            fObj.close()
            

        return arr


class ChollaSkeweryGlobal:
    '''
    Cholla Skewer y Global object
        Holds information regarding skewer along y direction including methods
            to grab data along skewer.
        For the y-skewer, in order to match up with skewers outputs from Cholla,
            we need to grab data in order of (j, i, k). For details of 
            implementation, see src/analysis/lya_statistics.cpp Line 1264

        Initialized with:
        - ChollaSkewerGlobalHead (ChollaSkewerGlobalHead): global head object
        - SnapPath (str): path to snapshot directory
        - ChollaGrid (ChollaGrid): grid object
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, ChollaGrid):
        self.skewGlobalHead = ChollaSkewerGlobalHead
        self.SnapPath = SnapPath

        self.ny_local = int(self.skewGlobalHead.n_los / self.skewGlobalHead.nlos_proc)

        # save all of the ChollaBoxHead objects that this skewer passes through
        boxheads = self.skewGlobalHead.nlos_proc * [None]

        # grab j, k coordinates
        j = self.skewGlobalHead.get_globalj()
        k = self.skewGlobalHead.get_globalk()

        for curr_nlosproc in range(self.skewGlobalHead.nlos_proc):
            # starting x coordinate for skewer within current LOS process
            ystart = int(self.ny_local * curr_nlosproc)
            boxhead = ChollaGrid.get_BoxHead_ijk(j, ystart, k)
            boxheads[curr_nlosproc] = boxhead

        self.boxheads = tuple(boxheads)

    
    def get_hydrodata(self, key, dtype=np.float32):
        '''
        Return a specific hydro dataset along skewer

        Args:
            key (str): key to access data from hydro hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        arr = np.zeros(self.skewGlobalHead.n_los, dtype=dtype)
        j_local = self.skewGlobalHead.skewLocalFaceHead.localface_joffset
        k_local = self.skewGlobalHead.skewLocalFaceHead.localface_koffset

        for boxhead in self.boxheads:
            box = ChollaBox(self.SnapPath, boxhead)

            # check key with boxhead having 0 y-offset
            if (boxhead.offset[1] == 0):
                assert box.check_hydrokey(key)

            fObj = h5py.File(box.get_hydrofPath(), 'r')
            arr[boxhead.offset[1] : boxhead.offset[1] + boxhead.local_dims[1]] = fObj.get(key)[j_local, :, k_local]
            fObj.close()


        return arr


class ChollaSkewerzGlobal:
    '''
    Cholla Skewer z Global object
        Holds information regarding skewer along z direction including methods
            to grab data along skewer.

        Initialized with:
        - ChollaSkewerGlobalHead (ChollaSkewerGlobalHead): global head object
        - SnapPath (str): path to snapshot directory
        - ChollaGrid (ChollaGrid): grid object
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, ChollaGrid):
        self.skewGlobalHead = ChollaSkewerGlobalHead
        self.SnapPath = SnapPath

        self.nz_local = int(self.skewGlobalHead.n_los / self.skewGlobalHead.nlos_proc)

        # save all of the ChollaBoxHead objects that this skewer passes through
        boxheads = self.skewGlobalHead.nlos_proc * [None]

        # grab j, k coordinates
        j = self.skewGlobalHead.get_globalj()
        k = self.skewGlobalHead.get_globalk()

        for curr_nlosproc in range(self.skewGlobalHead.nlos_proc):
            # starting x coordinate for skewer within current LOS process
            zstart = int(self.nz_local * curr_nlosproc)
            boxhead = ChollaGrid.get_BoxHead_ijk(j, k, zstart)
            boxheads[curr_nlosproc] = boxhead

        self.boxheads = tuple(boxheads)


    def get_hydrodata(self, key, dtype=np.float32):
        '''
        Return a specific hydro dataset along skewer

        Args:
            key (str): key to access data from hydro hdf5 file
            dtype (np type): (optional) numpy precision to use
        Returns:
            arr (arr): requested dataset
        '''

        arr = np.zeros(self.skewGlobalHead.n_los, dtype=dtype)
        j_local = self.skewGlobalHead.skewLocalFaceHead.localface_joffset
        k_local = self.skewGlobalHead.skewLocalFaceHead.localface_koffset

        for boxhead in self.boxheads:
            box = ChollaBox(self.SnapPath, boxhead)

            # check key with boxhead having 0 y-offset
            if (boxhead.offset[2] == 0):
                assert box.check_hydrokey(key)

            fObj = h5py.File(box.get_hydrofPath(), 'r')
            arr[boxhead.offset[2] : boxhead.offset[2] + boxhead.local_dims[2]] = fObj.get(key)[j_local, k_local, :]
            fObj.close()


        return arr




