import numpy as np
import os
from cholla_api.snap.ChollaSnap import ChollaSnap

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

    def make_movie(self, plot_fn, keys, imgfbase, imgftype='png', movie_nsnaps=None, valuecalcs=None, plt_kwargs=None):
        if movie_nsnaps is None:
            movie_nsnaps = range(self.totnSnap)
        if plt_kwargs is None:
            plt_kwargs = {}
        fnums = int(np.ceil(np.log10(len(movie_nsnaps))) + 1)
        for n, nsnap in enumerate(movie_nsnaps):
            ch_snap = self.createSnap(nsnap, keys=keys, load_data=True, snap_head=True)
            if valuecalcs is not None:
                ch_snap.calc_vals(valuecalcs)
            if plt_kwargs.get('save'):
                fnum = str(n).zfill(fnums)
                imgfout = f"{self.imgsPath}/{imgfbase}_{fnum}.{imgftype}"
                plt_kwargs['imgfout'] = imgfout
            plot_fn(ch_snap.data, ch_snap.head, plt_kwargs)

    def compare_plot(self, plot_fn, keys, nsnap1, nsnap2, imgfbase, imgftype='png', valuecalcs=None, plt_kwargs=None):
        snap1 = ChollaSnap(nsnap1, keys=keys, load_data=True, snap_head=False)
        snap2 = ChollaSnap(nsnap2, keys=keys, load_data=True, snap_head=False)
        if valuecalcs is not None:
            snap1.calc_vals(valuecalcs)
            snap2.calc_vals(valuecalcs)
        if plt_kwargs.get('save'):
            imgfout = f"{self.imgsPath}/{imgfbase}_{nsnap1:.0d}_{nsnap2:.0d}.{imgftype}"
            plt_kwargs['imgfout'] = imgfout
        plot_fn(snap1.data, snap1.head, snap2.data, snap2.head, plt_kwargs)


