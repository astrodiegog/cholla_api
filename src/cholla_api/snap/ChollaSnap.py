import numpy as np
import h5py


'''
My code has a couple of if(cholladatatype) where I do the same
    action depending on whether its gravity, particles, or hydro
    
    I could keep doing
        general function
        specific function 1, 2, 3
    
    I could also create an index / id for each data type loaded as var, not fixed.
        this way, I'm not passing strings but instead chsnap.datagrav_id
'''


class ChollaSnap:
    '''
    Load in the data onto a dictionary. Optionally also save the header
    '''

    def __init__(self, nSnap, dataDir, namebase, nBoxes, hydro=True, gravy=False, parts=False):
        self.dataDir = dataDir + "/" + str(nSnap)
        self.nSnap = nSnap
        self.namebase = namebase
        self.nBoxes = nBoxes
        self.head = None
        self.dims = None
        self.dims_loc = None
        
        self.data_hydro, self.data_gravy, self.data_parts = None, None, None
        
        if hydro:
            self.data_hydro = {}
        if gravy:
            self.data_gravy = {}
        if parts:
            self.data_parts = {}

    def key_datacheck(self, dataftype, key, key_str=None, raise_keyerr=False):
        '''
        check whether a key is inside the data
        '''
        cholla_data = None
        if dataftype == "hydro":
            cholla_data = self.data_hydro
        elif dataftype == "parts":
            cholla_data = self.data_parts
        elif dataftype == "gravy":
            cholla_data = self.data_gravy
        else:
            # handle error where file is alive or ftype not specified
            return -1
        
        if key_str is None:
            key_str = key
        
        try:
            cholla_data[key]
        except KeyError:
            err_message = '-- Error --\n'
            err_message += f'Unable to find {key_str} \n'
            err_message += f'Missing {key} in snapshot data \n'
            print(err_message)
            if raise_keyerr:
                raise
    
    def key_hydrodatacheck(self, key, key_str=None, raise_keyerr=False):
        self.key_datacheck("hydro", key, key_str=key_str, raise_keyerr=raise_keyerr)
    
    def key_partsdatacheck(self, key, key_str=None, raise_keyerr=False):
        self.key_datacheck("parts", key, key_str=key_str, raise_keyerr=raise_keyerr)
        
    def key_gravydatacheck(self, key, key_str=None, raise_keyerr=False):
        self.key_datacheck("gravy", key, key_str=key_str, raise_keyerr=raise_keyerr)
    
    def set_head(self):
        """
        set the header with the hydro file
        """
        nbox = 0
        fName = '{0}.{1}.{2}'.format(self.nSnap, self.namebase, nbox)
        fPath = self.dataDir + '/' + fName
        fObj_hydro = h5py.File(fPath, 'r')
        self.head = dict(fObj_hydro.attrs)
        fObj_hydro.close()
        self.dims = self.head['dims']
        
        # where there are more cells than boxes
        cells_over_nbox_mask = (self.dims / self.nBoxes)>1 
        test_modulos = self.dims[cells_over_nbox_mask]
        # if cells divide evenly into nBoxes, can save dims_loc as attr
        always_set_locals = False
        for test_modulo in test_modulos:
            always_set_locals = (test_modulo % self.nBoxes)>0

        if not always_set_locals:
            self.dims_loc = self.head['dims_local']
    
    ### get file path functions
    def get_hydro_fPath(self, nbox):
        fName = '{0}.{1}.{2}'.format(self.nSnap, self.namebase, nbox)
        return self.dataDir + '/' + fName
    
    def get_parts_fPath(self, nbox):
        fName = '{0}_particles.{1}.{2}'.format(self.nSnap, self.namebase, nbox)
        return self.dataDir + '/' + fName
    
    def get_gravy_fPath(self, nbox):
        fName = '{0}_gravity.{1}.{2}'.format(self.nSnap, self.namebase, nbox)
        return self.dataDir + '/' + fName
    
        
    ### loading data functions
    def load_hydro(self, keys):
        """
        load hydro data.
            use the offset and local dimensions to globally combine them
        """
        if self.data_hydro is None:
            return -1
        self.hydro_keys = keys
        for nbox in range(self.nBoxes):
            fPath = self.get_hydro_fPath(nbox)
            fObj = h5py.File(fPath, 'r')
            
            if self.dims_loc is None:
                curr_dimsloc = fObj.attrs['dims_local']
            else:
                curr_dimsloc = self.dims_loc
            curr_offset = fObj.attrs['offset']
            
            keyStart_x, keyStart_y, keyStart_z = curr_offset
            keyEnd_x, keyEnd_y, keyEnd_z = curr_offset + curr_dimsloc
            for key in self.hydro_keys:
                if nbox == 0:
                    self.data_hydro[key] = np.zeros(self.dims)
                curr_keydata = np.array(fObj[key]).reshape(curr_dimsloc)
                self.data_hydro[key][keyStart_x:keyEnd_x, keyStart_y:keyEnd_y, keyStart_z:keyEnd_z] = curr_keydata

            fObj.close()
            
    def load_parts(self, keys):
        """
        load particle data.
            we don't know the number of particles a priori, so we have to load in the data
            then create a np array of the correct size and index the values in
        """
        if self.data_parts is None:
            return -1
        
        self.parts_keys = keys
        density_bool = False
        density_key = "density"
        if density_key in keys:
            density_bool = True
            keys.remove(density_key)
        
        particle_data_tmp = {}
        tot_num_parts = 0
        # arrays that hold cumulative number of particles indexed by nbox
        keyStarts = np.zeros(self.nBoxes, dtype=int)
        keyEnds = np.zeros(self.nBoxes, dtype=int)
        
        for nbox in range(self.nBoxes):
            fPath = self.get_parts_fPath(nbox)
            fObj = h5py.File(fPath, 'r')
            
            for k, key in enumerate(keys):
                # keys can only be 1D arrays of posn and vels
                if nbox == 0:
                    # initialize a dynamic array (list) for each key
                    particle_data_tmp[key] = []
                
                ########
                # may have to change from local to global coord system?
                curr_keydata = np.array(fObj[key])
                ########
                
                particle_data_tmp[key].append(curr_keydata)
                
                if k == 0:
                    keyStarts[nbox] = tot_num_parts
                    tot_num_parts += curr_keydata.size
                    keyEnds[nbox] = tot_num_parts
            
            if (density_bool):
                # density is same as hydro
                if self.dims_loc is None:
                    curr_dimsloc = fObj.attrs['dims_local']
                else:
                    curr_dimsloc = self.dims_loc
                curr_offset = fObj.attrs['offset']

                keyStart_x, keyStart_y, keyStart_z = curr_offset
                keyEnd_x, keyEnd_y, keyEnd_z = curr_offset + curr_dimsloc
                if nbox == 0:
                    self.data_parts[density_key] = np.zeros(self.dims)
                curr_keydata = np.array(fObj[density_key]).reshape(curr_dimsloc)
                self.data_parts[density_key][keyStart_x:keyEnd_x, keyStart_y:keyEnd_y, keyStart_z:keyEnd_z] = curr_keydata

            fObj.close()
        
        
        for key in keys:
            ### combine onto one array for each key
            self.data_parts[key] = np.zeros(tot_num_parts)
            
            for nbox in range(self.nBoxes):
                keyStart, keyEnd = keyStarts[nbox], keyEnds[nbox]
                self.data_parts[key][keyStart:keyEnd] = particle_data_tmp[key][nbox]

    def load_gravy(self, keys):
        """
        load gravity data.
            very similar to particle loading but for gravity
            current only key for gravity is 'potential', but I leave it generalized
            in case future keys are added on to gravity
        """
        if self.data_gravy is None:
            return -1
        self.gravy_keys = keys
        
        gravity_data_tmp = {}
        tot_num_gravs = 0
        # arrays that hold cumulative number of particles indexed by nbox
        keyStarts = np.zeros(self.nBoxes, dtype=int)
        keyEnds = np.zeros(self.nBoxes, dtype=int)
        
        for nbox in range(self.nBoxes):
            fPath = self.get_gravy_fPath(nbox)
            fObj = h5py.File(fPath, 'r')
            
            for k, key in enumerate(keys):
                # keys can only be 1D arrays of posn and vels
                if nbox == 0:
                    # initialize a dynamic array (list) for each key
                    gravity_data_tmp[key] = []
                
                ########
                # may have to change from local to global coord system?
                curr_keydata = np.array(fObj[key])
                ########
                
                gravity_data_tmp[key].append(curr_keydata)
                
                if k == 0:
                    keyStarts[nbox] = tot_num_gravs
                    tot_num_gravs += curr_keydata.size
                    keyEnds[nbox] = tot_num_gravs

            fObj.close()
            
        for key in keys:
            ### combine onto one array for each key
            self.data_gravy[key] = np.zeros(tot_num_gravs)
            
            for nbox in range(self.nBoxes):
                keyStart, keyEnd = keyStarts[nbox], keyEnds[nbox]
                self.data_gravy[key][keyStart:keyEnd] = gravity_data_tmp[key][nbox]
    
    ### concatenate files
    def concat_files(self, concat_ftype):
        cholla_type = None
        fName = None
        if concat_ftype == "hydro":
            cholla_type = self.data_hydro
            fName = '{0}.{1}'.format(self.nSnap, self.namebase)
        elif concat_ftype == "parts":
            cholla_type = self.data_parts
            fName = '{0}_particles.{1}'.format(self.nSnap, self.namebase)
        elif concat_ftype == "gravy":
            cholla_type = self.data_gravy
            fName = '{0}_gravity.{1}'.format(self.nSnap, self.namebase)
        else:
            # handle error where file is alive or ftype not specified
            return -1
        
        fPath = self.dataDir + '/' + fName
        
        # need to handle OSError Unable to create file (file exists) case
        with h5py.File(fPath,'x') as fObj:
            
            for i, (key, val) in enumerate(self.head.items()):
                fObj.attrs.create(name=key, data=val, shape=val.shape, dtype=val.dtype)
            
            for i, (key, val) in enumerate(cholla_type.items()):
                fObj.create_dataset(name=key, shape=val.shape, dtype=val.dtype,
                                    data=val)
        
    def concat_hydro(self):
        self.concat_files("hydro")
        
    def concat_parts(self):
        self.concat_files("parts")
        
    def concat_gravy(self):
        self.concat_files("hydro")
        
    
    def calc_vals(self, cholla_calcs):
        """
        cholla_calcs is a list of ChollaValueCalc instances
        """
        for k, cholla_calc in enumerate(cholla_calcs):
            try:
                cholla_calc.set_reqs(self.data_hydro, self.head)
            except:
                print('Unable to set the key requirements for', cholla_calc.calc_fn.__name__)

            req_keys = cholla_calc.reqs
            try:
                for key in req_keys:
                    self.data_hydro[key]

                if cholla_calc.kwarg_fn is not None:
                    cholla_calc.set_kwargs(self.data_hydro, self.head)
            except KeyError:
                print('Unable to run function', cholla_calc.calc_fn.__name__)
                print('Missing a value from ', cholla_calc.reqs)
            else:
                cholla_calc.calculate(self.data_hydro)
                

class ChollaValueCalc:
    '''
    class that makes it easy to determine what is required to calculate that value
    '''

    def __init__(self, calc_fn, req_fn, kwarg_fn):
        self.calc_fn = calc_fn
        self.kwarg_fn = kwarg_fn
        self.req_fn = req_fn
        self.kwargs = {}

    def set_reqs(self, data, head):
        self.reqs = self.req_fn(data, head)

    def set_kwargs(self, data, head):
        self.kwarg_fn(self.kwargs, data, head)

    def calculate(self, data):
        self.calc_fn(data, **self.kwargs)


