import argparse
import os

import numpy as np
import h5py

from cholla_api.OTFanalysis.ChollaOnTheFlySkewers import ChollaOnTheFlySkewers
from cholla_api.OTFanalysis.ChollaOnTheFlyFluxPowerSpectrum import ChollaOnTheFlyFluxPowerSpectrum


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
        description="Compute and append power spectra")

    parser.add_argument("skewfname", help='Cholla skewer output file name', type=str)

    parser.add_argument("dlogk", help='Differential log of step-size for power spectrum k-modes', type=float)

    parser.add_argument('-c', '--combine', help='Whether to combine power spectrum from each axis',
                        action='store_true')

    parser.add_argument('-v', '--verbose', help='Print info along the way',
                        action='store_true')

    return parser



def P_k_calc(OTFSkewers, dlogk, combine=True, verbose=False, precision=np.float64):
    '''
    Calculate the mean transmitted flux power spectrum along each axis and save
        onto skewer output file

    Args:
        OTFSkewers (ChollaOnTheFlySkewers): skewers object, interacts with files
        dlogk (float): differential step in log k-space
        combine (bool): (optional) whether to combine power spectrum from each axis
        verbose (bool): (optional) whether to print important information
        precision (np type): (optional) numpy precision to use in calculations
    Returns:
        ...
    '''

    OTFFluxPowerSpectrum = ChollaOnTheFlyFluxPowerSpectrum(dlogk, OTFSkewers)

    # calculate power spectra in each direction
    k_x, P_k_x = OTFFluxPowerSpectrum.get_FPS_x(precision)
    k_y, P_k_y = OTFFluxPowerSpectrum.get_FPS_y(precision)
    k_z, P_k_z = OTFFluxPowerSpectrum.get_FPS_z(precision)

    # open file and append each power spectrum as new "PowerSpectrum" group
    with h5py.File(OTFSkewers.OTFSkewersfPath, 'r+') as fObj:
        PS_group_key = 'PowerSpectrum'
        if PS_group_key not in fObj.keys():
            if verbose:
                print(f'\t...initializing power spectrum group for file {OTFSkewers.OTFSkewersfPath}')

            fObj.create_group(PS_group_key)

        PS_x_avg_key = 'P_x(k_x)'
        k_x_edges_key = 'k_x_edges'
        PS_y_avg_key = 'P_y(k_y)'
        k_y_edges_key = 'k_y_edges'
        PS_z_avg_key = 'P_z(k_z)'
        k_z_edges_key = 'k_z_edges'

        # ensure each key is in the PowerSpectrum Group
        if PS_x_avg_key not in fObj[PS_group_key].keys():
            if verbose:
                print(f'\t...initializing empty power spectrum and k mode arrays in x-axis for file {OTFSkewers.OTFSkewersfPath}')

            PS_x_empty = np.zeros(P_k_x.shape, dtype=P_k_x.dtype)
            k_x_empty = np.zeros(k_x.shape, dtype=k_x.dtype)
            fObj[PS_group_key].create_dataset(k_x_edges_key, data=k_x_empty)
            fObj[PS_group_key].create_dataset(PS_x_avg_key, data=PS_x_empty)

        if PS_y_avg_key not in fObj[PS_group_key].keys():
            if verbose:
                print(f'\t...initializing empty power spectrum and k mode arrays in y-axis for file {OTFSkewers.OTFSkewersfPath}')

            PS_y_empty = np.zeros(P_k_y.shape, dtype=P_k_y.dtype)
            k_y_empty = np.zeros(k_y.shape, dtype=k_y.dtype)
            fObj[PS_group_key].create_dataset(k_y_edges_key, data=k_y_empty)
            fObj[PS_group_key].create_dataset(PS_y_avg_key, data=PS_y_empty)

        if PS_z_avg_key not in fObj[PS_group_key].keys():
            if verbose:
                print(f'\t...initializing empty power spectrum and k mode arrays in z-axis for file {OTFSkewers.OTFSkewersfPath}')
            PS_z_empty = np.zeros(P_k_z.shape, dtype=P_k_z.dtype)
            k_z_empty = np.zeros(k_z.shape, dtype=k_z.dtype)
            fObj[PS_group_key].create_dataset(k_z_edges_key, data=k_z_empty)
            fObj[PS_group_key].create_dataset(PS_z_avg_key, data=PS_z_empty)

        if verbose:
            print(f'\t...assigning power spectrum and k edges in each axis file {OTFSkewers.OTFSkewersfPath}')
        fObj[PS_group_key][k_x_edges_key][:] = k_x[:]
        fObj[PS_group_key][PS_x_avg_key][:] = P_k_x[:]

        fObj[PS_group_key][k_y_edges_key][:] = k_y[:]
        fObj[PS_group_key][PS_y_avg_key][:] = P_k_y[:]

        fObj[PS_group_key][k_z_edges_key][:] = k_z[:]
        fObj[PS_group_key][PS_z_avg_key][:] = P_k_z[:]

        if combine:
            # before combining power spectrum along each direction, make sure they're of same shape
            assert np.array_equal(P_k_x.shape, P_k_y.shape)
            assert np.array_equal(P_k_x.shape, P_k_z.shape)

            # also need to assert that k_xyz are all within some tolerance level
            # for now assume they're the same, and save k_x by default !

            # combined power spectrum
            P_k = (P_k_x + P_k_y + P_k_z) / 3.
            PS_avg_key = 'P(k)'
            k_edges_key = 'k_edges'
            if PS_avg_key not in fObj[PS_group_key].keys():
                if verbose:
                    print(f'\t...initializing empty power spectrum average and k mode arrays for file {OTFSkewers.OTFSkewersfPath}')
                PS_empty = np.zeros(P_k.shape, dtype=P_k.dtype)
                k_empty = np.zeros(k_x.shape, dtype=k_x.dtype)
                fObj[PS_group_key].create_dataset(k_edges_key, data=k_empty)
                fObj[PS_group_key].create_dataset(PS_avg_key, data=PS_empty)

            if verbose:
                print(f'\t...assigning average power spectrum and k edges for file {OTFSkewers.OTFSkewersfPath}')
            fObj[PS_group_key][k_edges_key][:] = k_x[:]
            fObj[PS_group_key][PS_avg_key][:] = P_k[:]

    return



def main():
    '''
    Compute the power spectrum and append to skewer file
    '''

    # Create parser
    parser = create_parser()

    # Save args
    args = parser.parse_args()

    if args.verbose:
        print("we're verbose in this mf !")
        print(f"--- We are looking at skewer file : {args.skewfname} ---")
        print(f"--- We are placing power spectra in dlogk : {args.dlogk} ---")
        if args.combine:
            print(f"--- We are saving power spectra in each direction (!) ---")
        else:
            print(f"--- We are NOT saving power spectra in each direction (!) ---")

    precision = np.float64

    # ensure that local optical depth is a dataset
    fObj = h5py.File(args.skewfname, 'r')
    local_opticaldepth_key = 'taucalc_local'
    assert local_opticaldepth_key in fObj['skewers_x'].keys()
    assert local_opticaldepth_key in fObj['skewers_y'].keys()
    assert local_opticaldepth_key in fObj['skewers_z'].keys()
    fObj.close()

    # ensure dlogk is reasonable
    assert args.dlogk > 0

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

    # calculate the power spectra
    P_k_calc(OTFSkewers, args.dlogk, args.combine, args.verbose, precision=np.float64)




if __name__=="__main__":
    main()


