import numpy as np
import os
import pathlib

from cholla_api.snap.ChollaSnap import ChollaSnap
from cholla_api.viz.ChollaViz import ChollaViz



class ChollaRun:

    def __init__(self, basePath, tot_nSnap, namebase, nBoxes, data_dir='/data', img_dir='/imgs', test_name=""):
        self.basePath = basePath
        self.dataPath = self.basePath + data_dir
        self.imgsPath = self.basePath + img_dir
        self.totnSnap = tot_nSnap
        self.namebase = namebase
        self.nBoxes = nBoxes
        self.test_name = test_name
        
        self.check_totnsnap()
        
    def check_totnsnap(self):
        num_datafiles = len(os.listdir(self.dataPath))
        
        if (num_datafiles != self.totnSnap * self.nBoxes):
            err_message = f'''
            The given number of snapshots ({self.totnSnap:.0f}) and number of boxes ({self.nBoxes:.0f})
            \t does not match the number of files in data directory ({num_datafiles:.0f})
            '''
            
            raise Exception(err_message)

    def createSnap(self, nSnap, keys=[], load_data=True, snap_head=False):
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
                
                
        
        
    def compare_plot(self, plot_fn, keys, nsnap1, nsnap2, imgfbase, imgftype='png', valuecalcs=None, plt_kwargs=None):
        '''
        need to re-design with ChollaViz module
        '''
        snap1 = ChollaSnap(nsnap1, keys=keys, load_data=True, snap_head=False)
        snap2 = ChollaSnap(nsnap2, keys=keys, load_data=True, snap_head=False)
        if valuecalcs is not None:
            snap1.calc_vals(valuecalcs)
            snap2.calc_vals(valuecalcs)
        if plt_kwargs.get('save'):
            imgfout = f"{self.imgsPath}/{imgfbase}_{nsnap1:.0d}_{nsnap2:.0d}.{imgftype}"
            plt_kwargs['imgfout'] = imgfout
        plot_fn(snap1.data, snap1.head, snap2.data, snap2.head, plt_kwargs)


