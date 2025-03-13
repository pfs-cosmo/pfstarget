'''

module for correcting for galactic extinction


'''
import os
import numpy as np 
from astropy.table import Table, join

# absorption coefficients of HSC filters
absorptionCoeff = {
    "g"    : 3.240,
    "r"    : 2.276,
    "i"    : 1.633,
    "z"    : 1.263,
    "y"    : 1.075,
    "n387" : 4.007,
    "n468" : 3.351,
    "n515" : 2.939,
    "n527" : 2.855,
    "n656" : 2.077,
    "n718" : 1.812,
    "n816" : 1.458,
    "n921" : 1.187,
    "n973" : 1.083,
    "n1010": 1.013,
    "i945" : 1.134,
}


def _extinction_correct(hsc, method='sfd98', release='s23b', zeropoint=True): 
    ''' apply correction for galactic extinction using different methods (SFD98,
    Zhou DESI) and zero-point photometry correction  


    comments: 
    * CHH (03/06/2025): We may want to separate the zero-point photometry correction from
        the galactic extinction for clarity. 
    '''
    if method == 'sfd98': 
        #for band in ['g', 'r', 'i', 'z', 'y']: 
        #    if 'a_{band}' not in hsc.keys(): 
        #       raise: ValueError("extinction `a_*` not included in photometry") 
        
        # E(B -V)
        a_g = hsc["a_g"]
        a_r = hsc["a_r"]
        a_i = hsc["a_i"]
        a_z = hsc["a_z"]
        a_y = hsc["a_y"]
        
        # correct magnitudes for dust extinction *and* zero-point offset
        g_mag = hsc["g_cmodel_mag"] - a_g
        r_mag = hsc["r_cmodel_mag"] - a_r
        i_mag = hsc["i_cmodel_mag"] - a_i
        z_mag = hsc["z_cmodel_mag"] - a_z
        y_mag = hsc["y_cmodel_mag"] - a_y
        
        if zeropoint: 
            # get photometry zero-point offsets from wide.stellar_sequence_offset
            # (see Issue #8 for details)  
            grizy_offset = _get_zeropoint_correct(hsc['tract'], hsc['patch'],
                                                  release=release) 

            g_mag -= grizy_offset[0]
            r_mag -= grizy_offset[1]
            i_mag -= grizy_offset[2]
            z_mag -= grizy_offset[3]
            y_mag -= grizy_offset[4]
        
    elif method == 'desi': 
        import healpy as hp
        
        nside = 512 # healpix nside hardcoded
        _hsc = Table() 
        _hsc['HPXPIXEL'] = hp.ang2pix(nside, hsc['ra'], hsc['dec'], lonlat=True, nest=False)

        # read dust model 
        desi_dust = Table.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                            'dat', 'desi_dust_gr_512.fits'))

        # get E(B-V) value based on healpixel  
        ebv_desi = join(_hsc, desi_dust['EBV_GR', 'EBV_SFD', 'HPXPIXEL'], join_type='left', keys='HPXPIXEL') 

        a_g = absorptionCoeff['g'] * ebv_desi['EBV_GR']
        a_r = absorptionCoeff['r'] * ebv_desi['EBV_GR']
        a_i = absorptionCoeff['i'] * ebv_desi['EBV_GR']
        a_z = absorptionCoeff['z'] * ebv_desi['EBV_GR']
        a_y = absorptionCoeff['y'] * ebv_desi['EBV_GR']

        # correct magnitudes for dust extinction
        g_mag = hsc["g_cmodel_mag"] - a_g 
        r_mag = hsc["r_cmodel_mag"] - a_r 
        i_mag = hsc["i_cmodel_mag"] - a_i 
        z_mag = hsc["z_cmodel_mag"] - a_z 
        y_mag = hsc["y_cmodel_mag"] - a_y 

        if zeropoint: # implemented based on Jingjing's suggestion in Issue #8

            # get photometry zero-point offsets from wide.stellar_sequence_offset
            grizy_offset = _get_zeropoint_correct(hsc['tract'], hsc['patch'],
                                                  release=release) 

            # calculate discrepancy between SFD and DESI 
            debv = (ebv_desi['EBV_SFD'] - ebv_desi['EBV_GR'])
            delta_a_g = absorptionCoeff['g'] * debv 
            delta_a_r = absorptionCoeff['r'] * debv
            delta_a_i = absorptionCoeff['i'] * debv
            delta_a_z = absorptionCoeff['z'] * debv
            delta_a_y = absorptionCoeff['y'] * debv

            # g_mag_offset (with DESI dust) = g_mag_offset + delta_a_g
            g_mag -= (grizy_offset[0] + delta_a_g)
            r_mag -= (grizy_offset[1] + delta_a_r)
            i_mag -= (grizy_offset[2] + delta_a_i)
            z_mag -= (grizy_offset[3] + delta_a_z)
            y_mag -= (grizy_offset[4] + delta_a_y)

    else:
        raise NotImplementedError

    return g_mag, r_mag, i_mag, z_mag, y_mag 


def _get_zeropoint_correct(tract, patch, release='s23b'): 
    ''' return g/r/i/z/y-band photometric zeropoint correction based on tract
    and patch. 
    '''
    if release != 's23b': 
        raise ValueError("zero-point correction only for S23B")
     
    # read pdr3_wide.stellar_sequence_offsets
    foffset = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                'dat', 's23b_stellar_offsets.csv.gz')
    offsets = Table.read(foffset, format='csv') 
        
    # match offsets to the input tracts and patches
    _tp = Table([tract, patch], names=['tract', 'patch'])
    mtable = join(_tp, offsets, join_type='left')

    output = np.array([np.array(mtable['g_mag_offset']), 
                       np.array(mtable['r_mag_offset']), 
                       np.array(mtable['i_mag_offset']), 
                       np.array(mtable['z_mag_offset']), 
                       np.array(mtable['y_mag_offset'])]) 
    return output 
