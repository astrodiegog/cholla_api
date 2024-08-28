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
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, chGrid, chMacroFlags):
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
            boxhead = chGrid.get_BoxHead_ijk(xstart, j, k)
            boxheads[curr_nlosproc] = boxhead
        
        self.boxheads = tuple(boxheads)
        self.MacroFlags = chMacroFlags


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
            box = ChollaBox(self.SnapPath, boxhead, self.MacroFlags)

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
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, chGrid, chMacroFlags):
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
            boxhead = chGrid.get_BoxHead_ijk(j, ystart, k)
            boxheads[curr_nlosproc] = boxhead

        self.boxheads = tuple(boxheads)
        self.MacroFlags = chMacroFlags

    
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
            box = ChollaBox(self.SnapPath, boxhead, self.MacroFlags)

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
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, ChollaSkewerGlobalHead, SnapPath, chGrid, chMacroFlags):
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
            boxhead = chGrid.get_BoxHead_ijk(j, k, zstart)
            boxheads[curr_nlosproc] = boxhead

        self.boxheads = tuple(boxheads)
        self.MacroFlags = chMacroFlags


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
            box = ChollaBox(self.SnapPath, boxhead, self.MacroFlags)

            # check key with boxhead having 0 y-offset
            if (boxhead.offset[2] == 0):
                assert box.check_hydrokey(key)

            fObj = h5py.File(box.get_hydrofPath(), 'r')
            arr[boxhead.offset[2] : boxhead.offset[2] + boxhead.local_dims[2]] = fObj.get(key)[j_local, k_local, :]
            fObj.close()


        return arr


class ChollaSkewerAnalysisHead:
    '''
    Cholla Skewer Analysis Head

    Holds information regarding a skewer analysis

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - nlos_global (int): number of line-of-sight global cells
        - nj_global (int): number of global cells along j-dimension
        - nk_global (int): number of global cells along k-dimension
        - nlos_proc (int): number of processes along line-of-sight
        - nj_proc (int): number of processes along j-dimension
        - nk_proc (int): number of processes along k-dimension
    '''
    def __init__(self, n_stride, nlos_global, nj_global, nk_global, nlos_proc, nj_proc, nk_proc):
        self.n_stride = n_stride
        self.nlos_global = nlos_global
        self.nlos_proc = nlos_proc
        self.nj_proc = nj_proc

        self.nj_local = int(nj_global / nj_proc)
        self.nk_local = int(nk_global / nk_proc)

        self.nSkewers_j = int(nj_global / self.n_stride)
        self.nSkewers_k = int(nk_global / self.n_stride)
        self.nSkewersTotal = int(self.nSkewers_j * self.nSkewers_k)

        self.nSkewerslocal_j = int(self.nj_local / self.n_stride)
        self.nSkewerslocal_k = int(self.nk_local / self.n_stride)
        self.nSkewersLocal = int(self.nSkewerslocal_j * self.nSkewerslocal_k)

        self.nFaces_j, self.nFaces_k = int(nj_proc), int(nk_proc)
        self.nFaces_tot = int(self.nFaces_j * self.nFaces_k)


    def get_facehead(self, global_id):
        '''
        Return the ChollaSkewerFaceHead object corresponding to skewer global id
            Faces are tiled first along k-axis, then j-axis

        Args:
            global_id (int): id of the global skewer
        Return:
            skewfacehead (ChollaSkewerFaceHead): FaceHead for this global skewer
        '''

        face_id = int(global_id // self.nSkewersLocal)
        face_joffset = int( (face_id % self.nFaces_j) * self.nj_local)
        face_koffset = int( (face_id // self.nFaces_j) * self.nk_local)

        skewfacehead = ChollaSkewerFaceHead(face_id, face_joffset, face_koffset)

        return skewfacehead

    
    def get_localfacehead(self, global_id):
        '''
        Return the ChollaSkewerLocalFaceHead object corresponding to skewer 
            global id
            Local skewers are tiled first along j-axis, then k-axis

        Args:
            global_id (int): id of the global skewer
        Return:
            skewlocalhead (ChollaSkewerLocalFaceHead): LocalHead for this global skewer
        '''

        local_id = int(global_id % self.nSkewersLocal)
        local_joffset = int( (local_id // self.nSkewerslocal_j) * self.n_stride)
        local_koffset = int( (local_id % self.nSkewerslocal_j) * self.n_stride)

        skewlocalhead = ChollaSkewerLocalFaceHead(local_id, local_joffset, 
                                                  local_koffset)

        return skewlocalhead

    def get_globalhead(self, global_id):
        '''
        Return the ChollaSkewerGlobalHead object corresponding to skewer global id

        Args:
            global_id (int): id of the global skewer
        Return:
            skewglobalhead (ChollaSkewerGlobalHead): GlobalHead for this global 
                skewer id
        '''

        facehead = self.get_facehead(global_id)
        localhead = self.get_localfacehead(global_id)

        skewglobalhead = ChollaSkewerGlobalHead(global_id, facehead, localhead, 
                                                self.nlos_global, self.nlos_proc)

        return skewglobalhead
       



class ChollaSkewerxAnalysis:
    '''
    Cholla Skewer x Analysis object
        Holds information regarding skewer analysis along x direction including
            methods to grab ChollaSkewerxGlobal objects

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - SnapPath (str): path to snapshot directory
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, n_stride, SnapPath, chGrid, chMacroFlags):
        self.skewHead = ChollaSkewerAnalysisHead(n_stride, 
                                                 chGrid.nx_global, chGrid.ny_global, chGrid.nz_global,
                                                 chGrid.nproc_x, chGrid.nproc_y, chGrid.nproc_z)
        self.SnapPath = SnapPath
        self.grid = chGrid
        self.MacroFlags = chMacroFlags

    def get_skewer(self, global_id):
        '''
        Returns Cholla Skewer x Global object corresponding to global skewer id

        Args:
            global_id (int): id of the global skewer
        Return:
            skewglobal (ChollaSkewerxGlobal): SkewerxGlobal for this global 
                skewer id
        '''

        skewxglobalhead = self.skewHead.get_globalhead(global_id)
        
        skewglobal = ChollaSkewerxGlobal(skewxglobalhead, self.SnapPath, self.grid, 
                                         self.MacroFlags)

        return skewglobal


class ChollaSkeweryAnalysis:
    '''
    Cholla Skewer y Analysis object
        Holds information regarding skewer analysis along y direction including
            methods to grab ChollaSkeweryGlobal objects
        For the y-skewers, in order to match up with skewers outputs from Cholla,
            we need to grab data in order of (j, i, k). For details of 
            implementation, see src/analysis/lya_statistics.cpp Line 1264

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - SnapPath (str): path to snapshot directory
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, n_stride, SnapPath, chGrid, chMacroFlags):
        self.skewHead = ChollaSkewerAnalysisHead(n_stride,
                                                 chGrid.ny_global, chGrid.nx_global, chGrid.nz_global,
                                                 chGrid.nproc_y, chGrid.nproc_x, chGrid.nproc_z)
        self.SnapPath = SnapPath
        self.grid = chGrid
        self.MacroFlags = chMacroFlags

    def get_skewer(self, global_id):
        '''
        Returns Cholla Skewer y Global object corresponding to global skewer id

        Args:
            global_id (int): id of the global skewer
        Return:
            skewglobal (ChollaSkeweryGlobal): SkeweryGlobal for this global 
                skewer id
        '''

        skewyglobalhead = self.skewHead.get_globalhead(global_id)

        skewglobal = ChollaSkeweryGlobal(skewyglobalhead, self.SnapPath, self.grid, 
                                         self.MacroFlags)

        return skewglobal


class ChollaSkewerzAnalysis:
    '''
    Cholla Skewer z Analysis object
        Holds information regarding skewer analysis along z direction including
            methods to grab ChollaSkeweryGlobal objects

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - SnapPath (str): path to snapshot directory
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''

    def __init__(self, n_stride, SnapPath, chGrid, chMacroFlags):
        self.skewHead = ChollaSkewerAnalysisHead(n_stride,
                                                 chGrid.nz_global, chGrid.nx_global, chGrid.ny_global,
                                                 chGrid.nproc_z, chGrid.nproc_x, chGrid.nproc_y)
        self.SnapPath = SnapPath
        self.grid = chGrid
        self.MacroFlags = chMacroFlags

    def get_skewer(self, global_id):
        '''
        Returns Cholla Skewer z Global object corresponding to global skewer id

        Args:
            global_id (int): id of the global skewer
        Return:
            skewglobal (ChollaSkewerzGlobal): SkewerzGlobal for this global
                skewer id
        '''

        skewzglobalhead = self.skewHead.get_globalhead(global_id)

        skewglobal = ChollaSkewerzGlobal(skewzglobalhead, self.SnapPath, self.grid,
                                         self.MacroFlags)

        return skewglobal


