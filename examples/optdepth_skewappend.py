import argparse
import os

import numpy as np
import h5py

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

    parser.add_argument("skewfname", help='Cholla skewer output file name', type=str)

    parser.add_argument('-l', '--local', help='Whether to store local optical depths',
                        action='store_true')

    parser.add_argument('-v', '--verbose', help='Print info along the way', 
                        action='store_true')

    return parser

def init_taucalc(OTFSkewers, verbose=False, local=False):
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
        verbose (bool): (optional) whether to print important information
        local (bool): (optional) whether to save local optical depths
    Returns:
        ...

    '''

    with h5py.File(OTFSkewers.OTFSkewersfPath, 'r+') as fObj:
        if verbose:
            print(f'\t...initializing optical depth calculations for file {OTFSkewers.OTFSkewersfPath}')

        OTFSkewers_lst = [OTFSkewers.get_skewersx_obj(),
                          OTFSkewers.get_skewersy_obj(),
                          OTFSkewers.get_skewersz_obj()]

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
            if ((local) and 'taucalc_local' not in fObj[skew_key].keys()):
                taucalc_local = np.zeros((OTFSkewers_i.OTFSkewersiHead.n_skews, OTFSkewers_i.OTFSkewersiHead.n_i), 
                                          dtype=np.float64)
                fObj[skew_key].create_dataset('taucalc_local', data=taucalc_local)


    if verbose:
        print("...initialization complete !")
    
    return


def taucalc(OTFSkewers_i, skewCosmoCalc, precision=np.float64, verbose=False, local=False):
    '''
    Calculate the effective optical depth for each skewer along an axis

    Args:
        OTFSkewers_i (ChollaOnTheFlySkewers_i): holds all skewer info along an axis
        skewCosmoCalc (ChollaSkewerCosmoCalculator): holds optical depth function
        precision (np type): (optional) numpy precision to use
        verbose (bool): (optional) whether to print important information
        local (bool): (optional) whether to save local optical depths
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
            if local:
                fObj[skew_key]['taucalc_local'][nSkewerID] = taus

            if ((verbose) and ( (fObj[skew_key].attrs['taucalc_prog'] // 0.1) > progress_tenperc) ):
                print(f"--- Completed {fObj[skew_key].attrs['taucalc_prog'] * 100 : .0f} % at skewer {nSkewerID:.0f} ---")
                progress_tenperc += 1

    if verbose:
        print("Effective optical depth calculation completed along ", OTFSkewers_i.OTFSkewersiHead.skew_key)

    return

def main():
    '''
    Append the array of median optical depths for a skewer file
    '''

    # Create parser
    parser = create_parser()

    # Save args
    args = parser.parse_args()

    if args.verbose:
        print("we're verbose in this mf !")
        print(f"--- We are looking at skewer file : {args.skewfname} ---")
        if args.local:
            print(f"--- We are saving local optical depths (!) ---")
        else:
            print(f"--- We are NOT saving local optical depths (!) ---")

    precision = np.float64
    
    _ = '''
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
    '''

    # convert relative path to skewer file name to absolute file path
    cwd = os.getcwd()
    if args.skewfname[0] != '/':
        relative_path = args.skewfname
        args.skewfname = cwd + '/' + relative_path

    # seperate the skewer output number and skewer directory
    skewfName = args.skewfname.split('/')[-1]
    nSkewerOutput = int(skewfName.split('_')[0])
    skewersdir = args.skewfname[:-(len(skewfName)+1)]

    # create ChollaOTFSkewers object
    OTFSkewers = ChollaOnTheFlySkewers(nSkewerOutput, skewersdir)

    # create cosmology header
    chCosmoHead = ChollaCosmologyHead(OTFSkewers.Omega_M, OTFSkewers.Omega_R, 
                                      OTFSkewers.Omega_K, OTFSkewers.Omega_L,
                                      OTFSkewers.w0, OTFSkewers.wa, OTFSkewers.H0)

    # create skew cosmo calc object
    snapHead = ChollaSnapHead(nSkewerOutput + 1) # snapshots are index-1
    snapHead.a = OTFSkewers.current_a

    # add progress attribute, boolean mask for whether tau is calculated, and tau itself
    init_taucalc(OTFSkewers, args.verbose, args.local)
        
    OTFSkewers_lst = [OTFSkewers.get_skewersx_obj(), OTFSkewers.get_skewersy_obj(),
                      OTFSkewers.get_skewersz_obj()]


    # complete calculation
    for i, OTFSkewers_i in enumerate(OTFSkewers_lst):
        if args.verbose:
            print(f"Starting calculation along axis {i:.0f}")

        if (i == 0):
            nlos = OTFSkewers.nx
            dx = OTFSkewers.dx
        elif (i == 1):
            nlos = OTFSkewers.ny
            dx = OTFSkewers.dy
        elif (i == 2):
            nlos = OTFSkewers.nz
            dx = OTFSkewers.dz

        skewCosmoCalc = ChollaSkewerCosmoCalculator(snapHead, chCosmoHead, nlos, dx, precision)
        taucalc(OTFSkewers_i, skewCosmoCalc, precision, args.verbose, args.local)


        

if __name__=="__main__":
    main()

