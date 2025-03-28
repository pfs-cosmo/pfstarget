import numpy as np 


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
