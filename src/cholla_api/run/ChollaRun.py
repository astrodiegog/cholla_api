import numpy as np
import os
import glob
import pathlib

from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.viz.ChollaViz import ChollaViz
from cholla_api.viz.ChollaVizCompare import ChollaVizCompare


class ChollaRun:
    '''
    Class that holds important information to manipulate and study a Cholla simulation run
    '''

    def __init__(self, basePath, namebase='h5', data_dir='/data', img_dir='/imgs', test_name=""):
        self.basePath = basePath
        self.dataPath = self.basePath + data_dir
        self.imgsPath = self.basePath + img_dir
        self.namebase = namebase
        self.totnSnap = len(glob.glob1(self.dataPath, f"*.{self.namebase}.0"))
        self.nBoxes = len(glob.glob1(self.dataPath, f"0.{self.namebase}.*"))
        self.test_name = test_name
        
        self.check_totnsnap()
        
    def check_totnsnap(self):
        '''
        raise exception if the total number of snapshots is incorrect
        '''
        num_datafiles = len(os.listdir(self.dataPath))
        
        if (num_datafiles != self.totnSnap * self.nBoxes):
            err_message = f'''
            The given number of snapshots ({self.totnSnap:.0f}) and number of boxes ({self.nBoxes:.0f})
            \t does not match the number of files in data directory ({num_datafiles:.0f})
            '''
            
            raise Exception(err_message)

    def createSnap(self, nSnap, keys=[], load_data=True, snap_head=False):
        '''
        creates a ChollaSnap instance
        params:
            nSnap (int): the snapshot number to load
            keys (list): list of keys to load from the dataset
            load_data (bool): whether to load the key data or not
            snap_head (bool): whether to keep the sim head with snapshot, or pass onto ChollaRun class
        '''
        if nSnap > self.totnSnap:
            print('Invalid snap number')
            return -1
        ch_snap = ChollaSnap(nSnap, self.dataPath, self.namebase, self.nBoxes)
        if load_data:
            ch_snap.load_data(keys)
        if not snap_head:
            self.head = dict(ch_snap.head)
            ch_snap.head = None
        return ch_snap

    def beg_vs_fin(self, keys, imgftype='png', test_name="", valuecalcs=None, plots_type=None, plt_kwargs=None):
        '''
        make plot comparing initial vs final conditions
        '''
        if plots_type is None:
            plots_type = ["density", "velocity", "pressure"]
        if plt_kwargs is None:
            plt_kwargs = {}

        ch_snap1 = self.createSnap(0, keys=keys, load_data=True, snap_head=True)
        ch_snap2 = self.createSnap(self.totnSnap-1, keys=keys, load_data=True, snap_head=True)

        if valuecalcs is not None:
            ch_snap1.calc_vals(valuecalcs)
            ch_snap2.calc_vals(valuecalcs)
        
        ch_comp = ChollaVizCompare(ch_snap1, ch_snap2, test_name=test_name, plt_kwargs=plt_kwargs)

        if ("density" in plots_type):
            if plt_kwargs.get("save"):
                imgfout = f"{self.imgsPath}/density_ic.{imgftype}"
                plt_kwargs["imgfout"] = imgfout
            ch_comp.density(plt_kwargs)
         
        if ("pressure" in plots_type):
            if plt_kwargs.get("save"):
                imgfout = f"{self.imgsPath}/pressure_ic.{imgftype}"
                plt_kwargs["imgfout"] = imgfout
            ch_comp.pressure(plt_kwargs)

        if ("velocity" in plots_type):
            if plt_kwargs.get("save"):
                imgfout = f"{self.imgsPath}/velocity_ic.{imgftype}"
                plt_kwargs["imgfout"] = imgfout
            ch_comp.velocity(plt_kwargs)



    def make_movie(self, keys, imgfbase, imgftype='png', test_name="", movie_nsnaps=None, valuecalcs=None, movie_plots=None, plt_kwargs=None):
        '''
        helper function that loops over each movie_nsnap and saves a figure
        '''
        if movie_plots is None:
            movie_plots = ["density", "velocity", "pressure"]
        if movie_nsnaps is None:
            movie_nsnaps = range(self.totnSnap)
        if plt_kwargs is None:
            plt_kwargs = {}
        fnums = int(np.ceil(np.log10(len(movie_nsnaps))) + 1)
        
        progress_arr = np.arange(1,10)*0.1
        if len(movie_nsnaps) < 10:
            # less than 10 nsnaps, just look at 30, 50, 80%
            progress_arr = np.array([0.3, 0.5, 0.8])
        progress_ind = np.array(progress_arr*len(movie_nsnaps), dtype=int)
        progress_ind_curr = 0
        progress_ind_final = progress_arr.size
            
        for n, nsnap in enumerate(movie_nsnaps):
            ch_snap = self.createSnap(nsnap, keys=keys, load_data=True, snap_head=True)
            if valuecalcs is not None:
                ch_snap.calc_vals(valuecalcs)
            if plt_kwargs.get('save'):
                fnum_str = str(n).zfill(fnums)
                img_fbase = f"{imgfbase}_{fnum_str}.{imgftype}"
                
            ch_viz = ChollaViz(ch_snap, test_name=test_name, plt_kwargs=plt_kwargs)
            
            if ("density" in movie_plots):
                if plt_kwargs.get('save'):
                    density_dir = f"{self.imgsPath}/density"
                    pathlib.Path(density_dir).mkdir(parents=True, exist_ok=True)
                    imgfout = density_dir + f"/{img_fbase}"
                    plt_kwargs['imgfout'] = imgfout
                ch_viz.density(plt_kwargs)
                
            if ("velocity" in movie_plots):
                if plt_kwargs.get('save'):
                    velocity_dir = f"{self.imgsPath}/velocity"
                    pathlib.Path(velocity_dir).mkdir(parents=True, exist_ok=True)
                    imgfout = velocity_dir + f"/{img_fbase}"
                    plt_kwargs['imgfout'] = imgfout
                ch_viz.velocity(plt_kwargs)
                
            if ("pressure" in movie_plots):
                if plt_kwargs.get('save'):
                    pressure_dir = f"{self.imgsPath}/pressure"
                    pathlib.Path(pressure_dir).mkdir(parents=True, exist_ok=True)
                    imgfout = pressure_dir + f"/{img_fbase}"
                    plt_kwargs['imgfout'] = imgfout
                ch_viz.pressure(plt_kwargs)
            
            if (n == progress_ind[progress_ind_curr]):
                curr_progress = progress_arr[progress_ind_curr]
                print(f"--- Progress: {curr_progress*100:.0f}% ---")
                
                progress_ind_curr += 1
                if progress_ind_curr >= progress_arr.size:
                    # done with progress bar stuff
                    progress_ind_curr = -1
                


