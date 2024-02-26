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

    def __init__(self, basePath, namebase='h5', data_dir='/data', img_dir='/imgs', test_name="", hydro=True, gravy=False, parts=False):
        self.basePath = basePath
        self.dataPath = self.basePath + data_dir
        self.imgsPath = self.basePath + img_dir
        self.namebase = namebase
        dirsin_data = next(os.walk(self.dataPath))[1] # pray no more files are in data_dir
        self.totnSnap = len(dirsin_data)
        self.nBoxes = len(glob.glob1(self.dataPath + "/0", f"0.{self.namebase}.*"))
        self.test_name = test_name
        
        self.hydro = hydro
        self.gravy = gravy
        self.parts = parts
        

    def createSnap(self, nSnap, keys={}, load_data=True, snap_head=False, hydro=True, gravy=False, parts=False):
        '''
        creates a ChollaSnap instance
        params:
            nSnap (int): the snapshot number to load
            keys (dict): dictionary holding list of keys to load from dataset
            load_data (bool): whether to load the key data or not
            snap_head (bool): whether to keep the sim head with snapshot, or pass onto ChollaRun class
        '''
        if nSnap > self.totnSnap:
            print('Invalid snap number')
            return -1
        ch_snap = ChollaSnap(nSnap, self.dataPath, self.namebase, self.nBoxes, hydro=hydro, gravy=gravy, parts=parts)
        ch_snap.set_head() # needed to load in data
            
        if load_data:
            if hydro:
                ch_snap.load_hydro(keys["hydro"])
            if gravy:
                ch_snap.load_gravy(keys["gravy"])
            if parts:
                ch_snap.load_parts(keys["parts"])
        
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

    

    def make_movie(self, keys, imgfbase="img", imgftype='png', movie_nsnaps=None, valuecalcs=None, plot_vals=None, plt_kwargs=None):
        '''
        helper function that loops over each movie_nsnap and saves a figure
        
        keys (dict): dictionary holding list of keys to load from dataset
        imgfbase (str) - prefix to the file name being saved
        imgftype (str) - file extension
        movie_nsnaps (arr) - list/arr of snapshot numbers to use
        valuecalcs (list) - list of chollavalue calcs to use
        plot_vals (list) - list of values that serve as keys into data structs
        plt_kwargs (dict) - dictionary of importatnt plotting key word args
        '''
        if plot_vals is None:
            plot_vals = ["density", "velocity", "pressure"]
        if movie_nsnaps is None:
            movie_nsnaps = range(self.totnSnap)
        if plt_kwargs is None:
            plt_kwargs = {}
        # number of digits allowed to specify snapshot number
        fnums = int(np.ceil(np.log10(len(movie_nsnaps))) + 1)
        
        # create a progress array tracker
        progress_arr = np.arange(1,10)*0.1
        if len(movie_nsnaps) < 10:
            # less than 10 nsnaps, just look at 30, 50, 80%
            progress_arr = np.array([0.3, 0.5, 0.8])
        progress_ind = np.array(progress_arr*len(movie_nsnaps), dtype=int)
        progress_ind_curr = 0
        progress_ind_final = progress_arr.size
            
        for n, nsnap in enumerate(movie_nsnaps):
            ch_snap = self.createSnap(nsnap, keys=keys, load_data=True, snap_head=True, hydro=self.hydro, gravy=self.gravy, parts=self.parts)
            if valuecalcs is not None:
                ch_snap.calc_vals(valuecalcs)
            if plt_kwargs.get('save'):
                fnum_str = str(n).zfill(fnums)
                img_fname = f"{imgfbase}_{fnum_str}.{imgftype}"
                
            ch_viz = ChollaViz(ch_snap, test_name=self.test_name, plt_kwargs=plt_kwargs)
            
            for plot_val in plot_vals:
                if plt_kwargs.get('save'):
                    curr_plt_dir = f"{self.imgsPath}/{plot_val}"
                    # make sure plotting directory is there
                    if not os.path.isdir(curr_plt_dir):
                        print(f"Creating directory {curr_plt_dir}")
                        pathlib.Path(curr_plt_dir).mkdir(parents=True, exist_ok=False)
                        
                    imgfout = curr_plt_dir + f"/{img_fname}"
                    plt_kwargs['imgfout'] = imgfout
                
                curr_plt_fn = ch_viz.plot_keys[plot_val]
                curr_plt_fn(plt_kwargs)
            
            if (n == progress_ind[progress_ind_curr]):
                curr_progress = progress_arr[progress_ind_curr]
                print(f"--- Progress: {curr_progress*100:.0f}% ---")
                
                progress_ind_curr += 1
                if progress_ind_curr >= progress_arr.size:
                    # done with progress bar stuff
                    progress_ind_curr = -1
                


