import numpy as np
import h5py


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


