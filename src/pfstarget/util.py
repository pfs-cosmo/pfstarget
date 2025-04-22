import os
import numpy as np 
from astropy.table import Table, join


def healpixelize(ra, dec, nside=128): 
    ''' given RA and Dec return healpix number count. This is to calculate
    target/random counts
    '''
    import healpy as hp 
    # total number of pixels
    npix = hp.nside2npix(nside)  
    
    hpix = hp.ang2pix(nside, 
                      np.radians(90.0 - dec), 
                      np.radians(ra))

    uhpix, nhpix = np.unique(hpix, return_counts=True)
    hp_map = np.zeros(npix)
    hp_map[uhpix] = nhpix

    return hp_map 


def patch_qa(tract, patch, release='s23b'): 
    if release != 's23b': 
        raise NotImplementedError("patch_qa only for S23B")
    
    # read pdr3_wide.patch_qa
    fpatch = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                'dat', 'patch_qa.s23b.csv.gz')
    patches = Table.read(fpatch, format='csv') 
    patches = patches[~(patches['tract'].mask | patches['patch'].mask)]

    # match patch to the input tracts and patches
    _tp = Table([tract, patch], names=['tract', 'patch'])
    mtable = join(_tp, patches, join_type='left')
    return mtable 
