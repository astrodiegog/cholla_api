import numpy as np
import h5py

class ChollaSnap:
    '''
    Load in the data onto a dictionary. Optionally also save the header
    '''

    def __init__(self, nSnap, dataDir, namebase, nBoxes):
        self.dataDir = dataDir
        self.nSnap = nSnap
        self.namebase = namebase
        self.nBoxes = nBoxes
        self.head = None
        self.dims = None
        self.dims_loc = None
        self.data = {}

    def key_datacheck(self, key, key_str=None, raise_keyerr=False):
        '''
        check whether a key is inside the data
        '''
        if key_str is None:
            key_str = key
        
        try:
            self.data[key]
        except KeyError:
            err_message = '-- Error --\n'
            err_message += f'Unable to find {key_str} \n'
            err_message += f'Missing {key} in snapshot data \n'
            print(err_message)
            if raise_keyerr:
                raise
        
    def load_data(self, keys, set_head=True):
        """
        load all data
        """
        self.keys = keys
        for nbox in range(self.nBoxes):
            fName = '{0}.{1}.{2}'.format(self.nSnap, self.namebase, nbox)
            fPath = self.dataDir + '/' + fName
            fObj = h5py.File(fPath, 'r')
            if nbox == 0:
                if set_head:
                    self.head = dict(fObj.attrs)
                self.dims = self.head['dims']
                
                # where there are more cells than boxes
                cells_over_nbox_mask = (self.dims / self.nBoxes)>1 
                test_modulos = self.dims[cells_over_nbox_mask]
                # if cells divide evenly into nBoxes, can save dims_loc ass attr
                always_set_locals = False
                for test_modulo in test_modulos:
                    always_set_locals = (test_modulo % self.nBoxes)>0

                if not always_set_locals:
                    self.dims_loc = self.head['dims_local']
                
            
            if self.dims_loc is None:
                curr_dimsloc = fObj.attrs['dims_local']
            else:
                curr_dimsloc = self.dims_loc
            curr_offset = fObj.attrs['offset']
            
            keyStart_x, keyStart_y, keyStart_z = curr_offset
            keyEnd_x, keyEnd_y, keyEnd_z = curr_offset + curr_dimsloc
            for key in self.keys:
                if nbox == 0:
                    self.data[key] = np.zeros(self.dims)
                curr_keydata = np.array(fObj[key]).reshape(curr_dimsloc)
                self.data[key][keyStart_x:keyEnd_x, keyStart_y:keyEnd_y, keyStart_z:keyEnd_z] = curr_keydata

            fObj.close()

    def calc_vals(self, cholla_calcs):
        """
        cholla_calcs is a list of ChollaValueCalc instances
        """
        for k, cholla_calc in enumerate(cholla_calcs):
            try:
                cholla_calc.set_reqs(self.data, self.head)
            except:
                print('Unable to set the key requirements for', cholla_calc.calc_fn.__name__)

            req_keys = cholla_calc.reqs
            try:
                for key in req_keys:
                    self.data[key]

                if cholla_calc.kwarg_fn is not None:
                    cholla_calc.set_kwargs(self.data, self.head)
            except KeyError:
                print('Unable to run function', cholla_calc.calc_fn.__name__)
                print('Missing a value from ', cholla_calc.reqs)
            else:
                cholla_calc.calculate(self.data)


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


