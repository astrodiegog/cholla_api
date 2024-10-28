import argparse
import os

import numpy as np
import h5py

from cholla_api.run.ChollaGrid import ChollaGrid
from cholla_api.analysis.ChollaCosmoCalculator import ChollaCosmologyHead
from cholla_api.analysis.ChollaSkewersCalc import ChollaSkewerCosmoCalculator
from cholla_api.snap.ChollaSnap import ChollaSnapHead
from cholla_api.OTFanalysis.ChollaOnTheFlySkewers import ChollaOnTheFlySkewers




###
# Create command line arg parser
###

def create_parser():
    '''
    Create a command line argument parser that grabs the number of nodes
        and the parameter text file. Allow for verbosity

    Args:
        ...
    Returns:
        ...
    '''

    parser = argparse.ArgumentParser(
        description="Compute and append optical depth")

    parser.add_argument("nodes", help='Number of MPI processes used in sim',
                        type=int)

    parser.add_argument("param", help='Cholla parameter text file', type=str)

    parser.add_argument('-v', '--verbose', help='Print info along the way', 
                        action='store_true')

    return parser

def init_taucalc(OTFSkewers, nStride, verbose=False):
    '''
    Initialize the calculation of the effective optical depth. For each skewers_i axis
        group, create three things:
        1. (attribute) progress for optical depth
        2. (dataset) boolean mask whether optical depth has been calculated for 
            a specific skewer
        3. (dataset) median of local optical depths for a skewer
    
    Current implementation assumes the same nStride along each axis

    Args:
        OTFSkewers (ChollaOnTheFlySkewers): holds OTF skewers specific info
        nStride (int): stride cell number between skewers
        verbose (bool): (optional) whether to print important information
    Returns:
        ...

    '''

    with h5py.File(OTFSkewers.OTFSkewersfPath, 'r+') as fObj:
        if verbose:
            print(f'\t...initializing optical depth calculations for file {OTFSkewers.OTFSkewersfPath}')

        OTFSkewers_lst = [OTFSkewers.get_skewersx_obj(nStride),
                          OTFSkewers.get_skewersy_obj(nStride),
                          OTFSkewers.get_skewersz_obj(nStride)]

        # add progress attribute, boolean mask for whether tau is calculated, and tau itself
        for i, OTFSkewers_i in enumerate(OTFSkewers_lst):
            if verbose:
                print(f"\t\t...initializing arrays and attributes along axis {i:.0f}")
            skew_key = OTFSkewers_i.OTFSkewersiHead.skew_key

            taucalc_bool = np.zeros(OTFSkewers_i.OTFSkewersiHead.n_skews, dtype=bool)
            taucalc = np.zeros(OTFSkewers_i.OTFSkewersiHead.n_skews, dtype=np.float64)

            if 'taucalc_prog' not in dict(fObj[skew_key].attrs).keys():
                fObj[skew_key].attrs['taucalc_prog'] = 0.
            if 'taucalc_bool' not in fObj[skew_key].keys():
                fObj[skew_key].create_dataset('taucalc_bool', data=taucalc_bool)
            if 'taucalc_eff' not in fObj[skew_key].keys():
                fObj[skew_key].create_dataset('taucalc_eff', data=taucalc)

    if verbose:
        print("...initialization complete !")
    
    return


def taucalc(OTFSkewers_i, skewCosmoCalc, precision=np.float64, verbose=False):
    '''
    Calculate the effective optical depth for each skewer along an axis

    Args:
        OTFSkewers_i (ChollaOnTheFlySkewers_i): holds all skewer info along an axis
        skewCosmoCalc (ChollaSkewerCosmoCalculator): holds optical depth function
        precision (np type): (optional) numpy precision to use
        verbose (bool): (optional) whether to print important information
    Returns:
        ...
    '''

    skew_key = OTFSkewers_i.OTFSkewersiHead.skew_key

    with h5py.File(OTFSkewers_i.fPath, 'r+') as fObj:
        curr_progress = fObj[skew_key].attrs['taucalc_prog']
        progress_tenperc = int(curr_progress // 0.1)
        if verbose:
            print(f"Starting calculations at {100 * curr_progress:.2f} % complete")

        # loop over each skewer
        for nSkewerID in range(OTFSkewers_i.OTFSkewersiHead.n_skews):
            # skip skewers whose optical depth already calculated
            if (fObj[skew_key]['taucalc_bool'][nSkewerID]):
                continue

            # grab skewer data & calculate effective optical depth
            OTFSkewer = OTFSkewers_i.get_skewer_obj(nSkewerID)
            vel = OTFSkewer.get_losvelocity(precision)
            densityHI = OTFSkewer.get_HIdensity(precision)
            temp = OTFSkewer.get_temperature(precision)
            taus = skewCosmoCalc.optical_depth_Hydrogen(densityHI, vel, temp, use_forloop=True)
            tau_eff = np.median(taus)

            # update attr, bool arr, and tau arr
            fObj[skew_key].attrs['taucalc_prog'] += (1. / OTFSkewers_i.OTFSkewersiHead.n_skews)
            fObj[skew_key]['taucalc_bool'][nSkewerID] = True
            fObj[skew_key]['taucalc_eff'][nSkewerID] = tau_eff

            if ((verbose) and ( (fObj[skew_key].attrs['taucalc_prog'] // 0.1) > progress_tenperc) ):
                print(f"--- Completed {fObj[skew_key].attrs['taucalc_prog'] * 100 : .0f} % at skewer {nSkewerID:.0f} ---")
                progress_tenperc += 1

    if verbose:
        print("Effective optical depth calculation completed along ", OTFSkewers_i.OTFSkewersiHead.skew_key)

    return

def main():
    '''
    Append the array of median optical depths for each skewer file
    '''

    # Create parser
    parser = create_parser()

    # Save args
    args = parser.parse_args()

    if args.verbose:
        print("we're verbose in this mf !")
        print(f"--- We are using {args.nodes:.0f} number of nodes ---")
        print(f"--- We are looking at param text file : {args.param} ---")

    # create set of all param values we are interested in saving
    params2grab = {'nx', 'ny', 'nz', 'xmin', 'ymin', 'zmin', 'xlen', 'ylen', 'zlen',
                    'H0', 'Omega_M', 'Omega_L', 'Omega_K', 'Omega_R', 'Omega_b', 'w0', 'wa',
                    'analysis_scale_outputs_file', 'skewersdir', 'lya_skewers_stride'}
    
    # read in params from param text file
    params = {}
    with open(args.param, 'r') as paramfile:
        for line in paramfile:
            # strip whitespace, then split by key-value pair assignment
            keyval_str = '='
            linesplit = line.strip().split(keyval_str)
            is_keyvalpair = len(linesplit) == 2
            if is_keyvalpair:
                key_str, val_str = linesplit
                if key_str in params2grab:
                    params[key_str] = val_str
    
    if len(params) != len(params2grab):
        print(f'--- MISSING FOLLOWING PARAMS IN PARAM TXT FILE {args.param} ---')
        for param in params2grab:
            if param not in params.keys():
                print('\t - ', param)
        print('--- PLEASE FIX... EXITING ---')
        exit()

    # convert relative to absolute paths
    paramsPaths = {'skewersdir', 'analysis_scale_outputs_file'}
    cwd = os.getcwd()
    for path_param in paramsPaths:
        if params[path_param][0] != '/':
            relative_path = params[path_param]
            params[path_param] = cwd + '/' + relative_path

    # find number of skewers
    with open(params['analysis_scale_outputs_file']) as scale_outputs:
        numSkews_scale = len(scale_outputs.readlines())
    numSkews_skewdir = len(os.listdir(params['skewersdir']))
    if numSkews_scale != numSkews_skewdir:
        print(f'--- NUMBER OF OUTPUT SKEWERS {numSkews_skewdir:.0f} DOESNT MATCH OUTPUT FILE {numSkews_scale} ---')
        print('--- PLEASE FIX... EXITING ---')
        exit()
    else:
        numSkewerOutputs = numSkews_skewdir

    if args.verbose:
        scaleOutputsPath = params['analysis_scale_outputs_file']
        skewDirPath = params['skewersdir']
        print(f'Peaked in {scaleOutputsPath} and {skewDirPath} to find we should expect {numSkewerOutputs:.0f} skewer outputs')

    # create ChollaGrid object
    grid = ChollaGrid(args.nodes, int(params['nx']), int(params['ny']), int(params['nz']),
                      float(params['xmin']), float(params['ymin']), float(params['zmin']),
                      float(params['xlen']) + float(params['xmin']),
                      float(params['ylen']) + float(params['ymin']),
                      float(params['zlen']) + float(params['zmin']))

    # assume cube grid, such that line-of-sight and cell width are the same, just use x info
    nlos, dx = int(params['nx']), grid.dx
    precision = np.float64

    # create cosmology header
    chCosmoHead = ChollaCosmologyHead(float(params['Omega_M']), float(params['Omega_R']),
                                      float(params['Omega_K']), float(params['Omega_L']),
                                      float(params['w0']), float(params['wa']), float(params['H0']))

    for nSkewerOutput in range(numSkewerOutputs):
        # create ChollaOTFSkewers object
        OTFSkewers = ChollaOnTheFlySkewers(nSkewerOutput, params['skewersdir'], grid)

        # create skew cosmo calc object
        snapHead = ChollaSnapHead(nSkewerOutput + 1) # snapshots are index-1
        snapHead.a = OTFSkewers.current_a
        skewCosmoCalc = ChollaSkewerCosmoCalculator(snapHead, chCosmoHead, nlos, dx, precision)

        init_taucalc(OTFSkewers, int(params['lya_skewers_stride']), args.verbose)
        
        OTFSkewers_lst = [OTFSkewers.get_skewersx_obj(int(params['lya_skewers_stride'])),
                          OTFSkewers.get_skewersy_obj(int(params['lya_skewers_stride'])),
                          OTFSkewers.get_skewersz_obj(int(params['lya_skewers_stride']))]


        # add progress attribute, boolean mask for whether tau is calculated, and tau itself
        for i, OTFSkewers_i in enumerate(OTFSkewers_lst):
            if args.verbose:
                print(f"Starting calculation along axis {i:.0f}")

            taucalc(OTFSkewers_i, skewCosmoCalc, precision, args.verbose)



        

if __name__=="__main__":
    main()

