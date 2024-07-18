'''
Cholla Grid
'''

import numpy as np
from cholla_api.data.ChollaBox import ChollaBoxHead


def greatest_prime_factor(num):
    '''
    Find the greatest prime factor of a number
    
    Args:
        num (int): number to factorize
    Returns:
        prime_factor (int): largest prime factor
    '''
    
    assert type(num) == int
    
    if ((num == 1) or (num == 2)):
        return num
    
    prime_factor = 2
    
    while (True):
        
        # keep dividing while evenly divisible
        while ((num % prime_factor) == 0):
            num /= prime_factor
        
        # cannot divide anymore
        if (num == 1):
            break
        
        # iterate to the next prime factor
        prime_factor += 1
    
    return prime_factor



class ChollaGrid:
    '''
    Cholla Grid object
        Holds global information about the domain and bounds of a simulation
            run. To complete domain decomposition, also need the number of
            processes simulation was ran with

        Initialized with:
        - nprocs (int): number of MPI processes
        - nx (int): number of cells in x dimension
        - ny (int): number of cells in y dimension
        - nz (int): number of cells in z dimension
        - xmin (float): minimum x coordinate
        - ymin (float): minimum y coordinate
        - zmin (float): minimum z coordinate
        - xmax (float): maximum x coordinate
        - ymax (float): maximum y coordinate
        - zmax (float): maximum z coordinate
    
    '''
    
    def __init__(self, nprocs, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax):

        if nprocs > 1:
            # odd number of processors not supported
            assert not nprocs % 2

        self.nprocs = nprocs

        # set number of cells globally
        self.nx_global, self.ny_global, self.nz_global = nx, ny, nz

        # set min, length, and size of cell in each dimension
        self.xmin = xmin
        self.Lx = xmax - self.xmin
        self.dx = self.Lx / self.nx_global

        self.ymin = ymin
        self.Ly = ymax - self.ymin
        self.dy = self.Ly / self.ny_global

        self.zmin = zmin
        self.Lz = zmax - self.zmin
        self.dz = self.Lz / self.nz_global

        # perform domain decomposition
        self.domain_decomp()

    
    def tile_decomposition(self):
        '''
        Tile the MPI processes in a block arrangement. Set the number
            of processes in each dimension as attributes

        Args:
            ...
        Returns:
            ...
        '''
        np_x, np_y, np_z = 1, 1, 1
        nproc_tmp = int(self.nprocs)

        # 1D Case - eZ
        if ((self.nz_global == 1) and (self.ny_global == 1)):
            np_x = nproc_tmp

            self.nproc_x = np_x
            self.nproc_y = np_y
            self.nproc_z = np_z
            return

        # get greatest prime factor of number of MPI processes
        n_gpf = greatest_prime_factor(nproc_tmp)

        # 2D Case - assign domain
        if (self.nz_global == 1):
            np_x = n_gpf
            np_y = nproc_tmp // np_x
    
            # make sure np_x > np_y
            if (np_x < np_y):
                temp = np_x
                np_x = np_y
                np_y = temp

            self.nproc_x = np_x
            self.nproc_y = np_y
            self.nproc_z = np_z
            return

        # 3D Case - while loop nprocs
        index = 0
        while (nproc_tmp > 1):
            n_gpf = greatest_prime_factor(nproc_tmp)
            nproc_tmp = nproc_tmp // n_gpf

            if ((index % 3) == 0):
                np_x = int(np_x * n_gpf)
            elif ((index % 3) == 1):
                np_y = int(np_y * n_gpf)
            else:
                np_z = int(np_z * n_gpf)

            index += 1

        # ensure nx > ny > nz order
        if (np_z > np_y):
            temp = np_y
            np_y = np_z
            np_z = temp
        if (np_y > np_x):
            temp = np_x
            np_x = np_y
            np_y = temp
        if (np_z > np_y):
            temp = np_y
            np_y = np_z
            np_z = temp
        
        self.nproc_x = np_x
        self.nproc_y = np_y
        self.nproc_z = np_z

        return


    def domain_decomp(self):
        '''
        Tile the processes, set local subdomain sizes, and set starting 
            coords as attributes
    
        WARNING: assumes even split of global bounds
        
        Args:
            ...
        Returns:
            ...
        '''

        self.tile_decomposition()
    
        # set local x
        n = self.nx_global % self.nproc_x
        if (not n):
            self.nx_local = self.nx_global // self.nproc_x
    
        # set local y
        n = self.ny_global % self.nproc_y
        if (not n):
            self.ny_local = self.ny_global // self.nproc_y
    
        # set local z
        n = self.nz_global % self.nproc_z
        if (not n):
            self.nz_local = self.nz_global // self.nproc_z
        
        # set indices for each process
        n = 0
        self.index_x = np.zeros(self.nprocs, dtype=int)
        self.index_y = np.zeros(self.nprocs, dtype=int)
        self.index_z = np.zeros(self.nprocs, dtype=int)
        
        for k in range(self.nproc_z):
            for j in range(self.nproc_y):
                for i in range(self.nproc_x):
                    self.index_x[n] = i
                    self.index_y[n] = j
                    self.index_z[n] = k
                    n += 1

        return

    def get_BoxHead(self, nBox):
        '''
        Return the ChollaBoxHead object of this process number

        Args:
            nBox (int): process ID to use
        Return:
            boxhead (ChollaBoxHead): BoxHead for this object
        '''

        # calculate number of cell offsets
        offset_x = self.index_x[nBox] * self.nx_local
        offset_y = self.index_y[nBox] * self.ny_local
        offset_z = self.index_z[nBox] * self.nz_local

        # create object
        boxhead = ChollaBoxHead(nBox, offset_x, offset_y, offset_z,
                                self.nx_local, self.ny_local,
                                self.nz_local)
        # assign coordinates
        boxhead.set_coords(self.xmin, self.ymin, self.zmin, self.dx,
                               self.dy, self.dz)

        return boxhead


    def get_BoxHeads(self):
        '''
        Return a list of BoxHead objects

        Args:
            ...
        Returns:
            (tuple): BoxHeads for each process
        '''

        boxheads = []

        for nproc in range(self.nprocs):
            boxhead = self.get_BoxHead(nproc)

            # append onto list
            boxheads.append(boxhead)

        return tuple(boxheads)
    


