"""
Unpublished method (publication in progress); to try this out,
make sure you have numpy, nibabel, and dipy.

python cci.py name_of_any_trackvis_file.trk

"""

import numpy as np
import nibabel as nib
from dipy.tracking.metrics import downsample
from dipy.tracking.distances import bundles_distances_mdf
from dipy.tracking.utils import length

def save_cci_trkfile(tvis_hdr, sls, score_mtrx, save_path):
    """ Saves the cluster confidence index (cci), which is an
    estimation of the support a set of streamlines gives to
    a particular pathway, to a trk file. This enables the user
    to filter streamlines by cci in Trackvis. 
    (This is like a density threshold applied to pathways instead of voxels)

    Ex: A single streamline with no others in the dataset
    following a similar pathway has a low cci. A streamline
    in a bundle of 100 streamlines that follow similar
    pathways has a high cci.

    See: insert hypothetical publication here
    Based on MDF distance from Garyfallidis et al. 2012, see Quickbundles

    Parameters
    ----------
    tvis_hdr : header from a trackvis file
    streamlines : list of 2D (N, 3) arrays
        A sequence of streamlines of length N (# streamlines)
    score_mtrx : 2D array (N,1)
        An array of cci values in order of streamlines
    save_path : string
        The path where a trk file containing cci values will be saved

    Returns
    -------
    Returns nothing: the streamlines are saved at the specified path

    References
    ----------
    .. [Garyfallidis12] Garyfallidis E. et al., QuickBundles a method for
                        tractography simplification, Frontiers in Neuroscience,
                        vol 6, no 175, 2012.

    """
    #note the cci in the header as a streamline property
    newhdr=tvis_hdr.copy()
    newhdr['n_properties']=1
    proplist=newhdr['property_name']
    proplist[0]='Cluster Confidence'
    newhdr['property_name']=proplist

    #save the streamlines with the cci as a property
    newtrk = ((sl, None, score_mtrx[i]) for i,sl in enumerate(sls))
    nib.trackvis.write(save_path,newtrk,newhdr)

def calculate_cci(streamlines, max_mdf=5, min_sl_length=20, subsample=12, pow=1):
    """ Computes the cluster confidence index (cci), which is an
    estimation of the support a set of streamlines gives to
    a particular pathway.

    Ex: A single streamline with no others in the dataset
    following a similar pathway has a low cci. A streamline
    in a bundle of 100 streamlines that follow similar
    pathways has a high cci.

    See: insert hypothetical publication here
    Based on MDF distance from Garyfallidis et al. 2012,

    Parameters
    ----------
    streamlines : list of 2D (N, 3) arrays
        A sequence of streamlines of length N (# streamlines)
    max_mdf : int
        The maximum MDF distance (mm) that will be considered a
        "supporting" streamline and included in cci calculation
    min_sl_length : int
        The minimum length (in mm) of streamline that will be
        considered as a "supporting" streamline. Streamlines
        shorter than min_sl_length will not be included in the
        final trk file
    subsample: int
        The number of points that are considered for each streamline
        in the calculation. To save on calculation time, each
        streamline is subsampled to subsampleN points.
    pow: int
        The power to which the MDF distance for each streamline
        will be raised to determine how much it contributes to
        the cci. High values of pow make the contribution value
        degrade much faster. Example: a streamline with 5mm MDF
        similarity contributes 1/5 to the cci if pow is 1, but
        only contributes 1/5^2 = 1/25 if pow is 2.

    Returns
    -------
    Returns nothing: the streamlines are saved at the specified path

    References
    ----------
    .. [Garyfallidis12] Garyfallidis E. et al., QuickBundles a method for
                        tractography simplification, Frontiers in Neuroscience,
                        vol 6, no 175, 2012.


    """

    #remove any streamlines that are shorter than the min_sl_length
    lengths = list(length(streamlines))
    long_streamlines = []
    for n,j in enumerate(streamlines):
        if len(j) > subsample: #this prevents an error... I think it didn't like evaluating length of very short streamlines
            if lengths[n]>min_sl_length:
                long_streamlines.append(j)

    #calculate the pairwise MDF distance between all streamlines in dataset that are longer than min_sl_length
    subsamp_sls = [downsample(sl, subsample) for sl in long_streamlines]

    cci_score_mtrx = np.zeros([len(subsamp_sls)])

    for i,sl in enumerate(subsamp_sls):
        mdf_mtrx = bundles_distances_mdf([subsamp_sls[i]], subsamp_sls)
        mdf_mtrx_oi = (mdf_mtrx > 0) & (mdf_mtrx < max_mdf) & ~ np.isnan(mdf_mtrx)
        mdf_mtrx_oi_only = mdf_mtrx[mdf_mtrx_oi]
        cci_score = np.sum(np.divide(1, np.power(mdf_mtrx_oi_only, pow)))
        cci_score_mtrx[i] = cci_score

    return cci_score_mtrx, long_streamlines

def calc_save_cci_trkfile(path_to_trkfile, save_path, max_mdf=5, min_sl_length=20, subsample=12, pow=1):
    """ Saves the cluster confidence index (cci), which is an
    estimation of the support a set of streamlines gives to
    a particular pathway, to a trk file. This enables the user
    to filter streamlines by cci in Trackvis.

    Ex: A single streamline with no others in the dataset
    following a similar pathway has a low cci. A streamline
    in a bundle of 100 streamlines that follow similar
    pathways has a high cci.

    See: insert hypothetical publication here
    Based on MDF distance from Garyfallidis et al. 2012, see Quickbundles

    Parameters
    ----------
    path_to_trkfile: string
        The path of a trk file you want to calculate cci on
    save_path : string
        The path where a trk file containing cci values will be saved

    Returns
    -------
    Returns nothing: the streamlines are saved at the specified path

    References
    ----------
    .. [Garyfallidis12] Garyfallidis E. et al., QuickBundles a method for
                        tractography simplification, Frontiers in Neuroscience,
                        vol 6, no 175, 2012.

    """

    trk, hdr = nib.trackvis.read(path_to_trkfile)
    sls = [item[0] for item in trk]
    cci_array, sls_oi = calculate_cci(sls, max_mdf, min_sl_length, subsample, pow)
    save_cci_trkfile(hdr, sls_oi, cci_array, save_path)

if __name__ == '__main__':
    import os
    import sys
    trkloc = os.path.abspath(sys.argv[1])
    putloc = trkloc.replace('.trk', '_cci.trk')
    calc_save_cci_trkfile(trkloc, putloc, min_sl_length=20, subsample=8)

    print "Now try opening the new trkfile in trackvis and filtering using Cluster Confidence on the menu"
