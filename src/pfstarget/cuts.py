'''

module for target selection 


'''
import numpy as np 

from . import extinction as E


def isCosmology(objects, star_galaxy_cut=-0.15, magnitude_cut=22.5,
                g_r_cut=0.15, color_slope=2.0, color_yint=-0.15):
    ''' Select targets for the PFS Cosmology Survey

    args: 
        objects: HSC objects 

    kwargs: 
        star_galaxy_cut: cut for the star-galaxy separation. (Default: -0.15) 

    
    references:
    ----------
    * DESI elg target selection https://github.com/desihub/desitarget/blob/main/py/desitarget/cuts.py#L646

    '''
    # mask
    is_mask = masking(objects) 

    # quality cuts
    is_quality = quality_cuts(objects)

    # star-galaxy separation
    is_galaxy = star_galaxy(objects, cut=star_galaxy_cut)

    # color cut
    is_color = color_cut(objects, magnitude_cut=magnitude_cut) 

    return ~is_mask & is_quality & is_galaxy & is_color 


def color_cut(objects, magnitude_cut=22.5, g_r_cut=0.15, color_slope=2.0, color_yint=-0.15): 
    ''' impose color cut to select ELG within the redshift range of 0.6 < z <
    2.4
    '''
    # faint magnitude cut 
    cuts = objects['G_MAG'] > magnitude_cut  
    
    # color cut 
    cuts &= (((objects['G_MAG'] - objects['R_MAG']) < 0.15) | # g-r cut (for 1.6 < z < 2.4 ELGs) 
             ((objects['I_MAG'] - objects['Z_MAG']) 
              > color_slope * (objects['G_MAG'] - objects['R_MAG'])- color_yint))
    return cuts 


def quality_cuts(objects): 
    ''' quality cuts. Currently we only impose a magnitude depend cut on g-band
    magnitude error.
    '''
    # has g, r, i, y-band magnitudes 
    cuts = (np.isfinite(objects['G_MAG']) & 
            np.isfinite(objects['R_MAG']) & 
            np.isfinite(objects['I_MAG']) & 
            np.isfinite(objects['Z_MAG']))

    # g-band magnitude error cut
    cuts &= (objects['G_ERR'] < objects['G_MAG'] * 0.05 - 1.1) 

    # cut on extreme colors (these values are preliminary --- revisit this in detail)
    cuts &= ((objects['G_MAG'] - objects['R_MAG'] > -1) & 
             (objects['I_MAG'] - objects['Z_MAG'] > -1))
    return cuts


def star_galaxy(objects, cut=-0.15): 
    ''' star-galaxy separation. Currently PSF star-galaxy separation is:

    [measured i-band cmodel magnitude] - [measured i-band PSF magnitude] < cut

    This cut essentially determines whether the object is extended and has
    significant light outside of the PSF profile. 
    '''
    return (objects['I_MEAS_CMODEL_MAG'] - objects['I_MEAS_PSF_MAG'] < cut)


def masking(objects): 
    ''' masks including mask around bright objects from ghost, halo, blooming

    see https://hscla.mtk.nao.ac.jp/doc/wp-content/uploads/2020/12/brightstarmask.pdf
    for additional details on the masking. 

    return: 
        boolean array that specifies the objects that are *within* the mask 
    '''
    _mask = objects['I_MASK_HALO'].astype(bool) # within bright star halo 
    _mask &= objects['I_MASK_GHOST'].astype(bool) # larger circular radius around the bright star 
    _mask &= objects['I_MASK_BLOOMING'].astype(bool) 
    return _mask 


def _prepare_hsc(hsc, dust_extinction='sfd98'): 
    ''' prepare HSC imaging data for target selection 

    args:
        hsc : object 
    
    return: 
        objects: structured numpy array of HSC objects with relevant columns
        for target selection 
    '''
    dtype = [('OBJID', '<i8'), 
             ('RA', 'f4'), 
             ('DEC', 'f4'), 
             ('G_MAG', 'f4'), 
             ('R_MAG', 'f4'), 
             ('I_MAG', 'f4'), 
             ('Z_MAG', 'f4'), 
             ('Y_MAG', 'f4'), 
             ('G_ERR', 'f4'), 
             ('R_ERR', 'f4'), 
             ('I_ERR', 'f4'), 
             ('Z_ERR', 'f4'), 
             ('Y_ERR', 'f4'), 
             ('I_MEAS_CMODEL_MAG', 'f4'), 
             ('I_MEAS_PSF_MAG', 'f4'), 
             ('I_MASK_HALO', 'i'), 
             ('I_MASK_GHOST', 'i'), 
             ('I_MASK_BLOOMING', 'i')
             ]

    objects = np.zeros(len(hsc), dtype=dtype)

    # grizy magnitudes corrections for galactic dust extinction 
    _g, _r, _i, _z, _y = E._extinction_correct(hsc, method='sfd'): 
    objects['G_MAG'] = _g
    objects['R_MAG'] = _r
    objects['I_MAG'] = _i
    objects['Z_MAG'] = _z
    objects['Y_MAG'] = _y

    # g-band CMODEL magnitude used for quality cuts 
    objects['G_ERR'] = hsc["g_cmodel_mag_err"]
    objects['R_ERR'] = hsc["r_cmodel_mag_err"]
    objects['I_ERR'] = hsc["i_cmodel_mag_err"]
    objects['Z_ERR'] = hsc["z_cmodel_mag_err"]
    objects['Y_ERR'] = hsc["y_cmodel_mag_err"]

    # i-band measured CMODEL and PSF magnitudes for star-galaxy separation 
    objects['I_MEAS_CMODEL_MAG'] = hsc["i_meas_cmodel_mag"]
    # TEMPORARILY INCLUDING TYPO FIX THIS BACK TO PSF LATER
    objects['I_MEAS_PSF_MAG'] = hsc["i_meas_pfs_mag"] 

    # i-band mask for bright stars
    objects['I_MASK_HALO']      = hsc["i_mask_brightstar_halo"].astype(int)
    objects['I_MASK_GHOST']     = hsc["i_mask_brightstar_ghost"].astype(int)
    objects['I_MASK_BLOOMING']  = hsc["i_mask_brightstar_blooming"].astype(int)

    # additional columns 
    objects['OBJID'] = hsc['object_id']  # object id 
    objects['RA'] = hsc['ra'] # RA from i-band measurement
    objects['DEC'] = hsc['dec'] # RA from i-band measurement

    return objects
