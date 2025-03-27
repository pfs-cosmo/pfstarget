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
    cuts = objects['I_MAG'] > magnitude_cut  # I_MAG <= 24 is hardcoded in the SQL 
    
    # color cut 
    cuts &= (((objects['G_MAG'] - objects['R_MAG']) < g_r_cut) | # g-r cut (for 1.6 < z < 2.4 ELGs) 
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

    # psf_flux_flag: psf_flux is required by the observatory
    cuts &= ((~objects['G_PSF_FLAG']) & 
             (~objects['R_PSF_FLAG']) & 
             (~objects['I_PSF_FLAG']) &
             (~objects['Z_PSF_FLAG']))

    # g-band magnitude error cut
    cuts &= (objects['G_ERR'] < objects['G_MAG'] * 0.05 - 1.1) 

    # not skipped by deblender (this is a very very small fraction of objects)
    cuts &= (~objects['DEBLEND_SKIPPED'])

    # low surface brightness object cut (ref. Li Xiangchong et al. 2022, Table
    # 2; also Issue #10) 
    cuts &= ((objects['I_APFLUX10_MAG'] <= 25.5) & 
             (~objects['I_APFLUX10_FLAG']))

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
    _mask = objects['i_mask_halo'].astype(bool) # within bright star halo 
    _mask &= objects['i_mask_ghost'].astype(bool) # larger circular radius around the bright star 
    _mask &= objects['i_mask_blooming'].astype(bool) 
    return _mask 


def _prepare_hsc(hsc, dust_extinction='sfd98', release='s23b', zeropoint=true): 
    ''' prepare hsc imaging data for target selection 

    args:
        hsc : object 
            astropy.table or some structured array with hsc data 

        dust_extinction : str
            string specifying the galactic dust extinction model 
            (default: 'sfd98) 

        release : str
            string specifying the hsc data release (default: s23b) 
    
    return: 
        objects: structured numpy array of hsc objects with relevant columns
        for target selection 
    '''
    dtype = [('objid', '<i8'), 
             ('ra', 'f4'), 
             ('dec', 'f4'), 
             ('g_mag', 'f4'), 
             ('r_mag', 'f4'), 
             ('i_mag', 'f4'), 
             ('z_mag', 'f4'), 
             ('y_mag', 'f4'), 
             ('g_err', 'f4'), 
             ('r_err', 'f4'), 
             ('i_err', 'f4'), 
             ('z_err', 'f4'), 
             ('y_err', 'f4'), 
             ('i_meas_cmodel_mag', 'f4'), 
             ('i_meas_psf_mag', 'f4'), 
             ('i_mask_halo', 'i'), 
             ('i_mask_ghost', 'i'), 
             ('i_mask_blooming', 'i'),
             ('g_psf_flag', bool), 
             ('r_psf_flag', bool),
             ('i_psf_flag', bool), 
             ('z_psf_flag', bool), 
             ('deblend_skipped', bool), 
             ('i_apflux10_mag', 'f4'),  
             ('i_apflux10_flag', bool)
             ]

    objects = np.zeros(len(hsc), dtype=dtype)

    # grizy magnitudes corrections for galactic dust extinction 
    _g, _r, _i, _z, _y = e._extinction_correct(hsc, method=dust_extinction, 
                                               release=release,
                                               zeropoint=zeropoint)
    objects['g_mag'] = _g
    objects['r_mag'] = _r
    objects['i_mag'] = _i
    objects['z_mag'] = _z
    objects['y_mag'] = _y

    # g-band cmodel magnitude used for quality cuts 
    objects['g_err'] = hsc["g_cmodel_mag_err"]
    objects['r_err'] = hsc["r_cmodel_mag_err"]
    objects['i_err'] = hsc["i_cmodel_mag_err"]
    objects['z_err'] = hsc["z_cmodel_mag_err"]
    objects['y_err'] = hsc["y_cmodel_mag_err"]

    # i-band measured cmodel and psf magnitudes for star-galaxy separation 
    objects['i_meas_cmodel_mag'] = hsc["i_meas_cmodel_mag"]
    # temporarily including typo fix this back to psf later
    objects['i_meas_psf_mag'] = hsc["i_meas_psf_mag"] 

    # i-band mask for bright stars
    objects['i_mask_halo']      = hsc["i_mask_brightstar_halo"].astype(int)
    objects['i_mask_ghost']     = hsc["i_mask_brightstar_ghost"].astype(int)
    objects['i_mask_blooming']  = hsc["i_mask_brightstar_blooming"].astype(int)
    
    # psf flag used for quality cut 
    objects['g_psf_flag'] = hsc['g_psf_flag']
    objects['r_psf_flag'] = hsc['r_psf_flag']
    objects['i_psf_flag'] = hsc['i_psf_flag']
    objects['z_psf_flag'] = hsc['z_psf_flag']
    
    # skipped by deblender 
    objects['deblend_skipped'] = hsc['deblend_skipped']

    # aperture flux 
    objects['i_apflux10_mag']    = hsc['i_apertureflux_10_mag']
    objects['i_apflux10_flag']   = hsc['i_apertureflux_10_flag']

    # additional columns 
    objects['objid'] = hsc['object_id']  # object id 
    objects['ra'] = hsc['ra'] # ra from i-band measurement
    objects['dec'] = hsc['dec'] # ra from i-band measurement
             
    return objects


def random_masking(randoms): 
    ''' masks including mask around bright objects from ghost, halo, blooming
    for random catalog 

    see https://hscla.mtk.nao.ac.jp/doc/wp-content/uploads/2020/12/brightstarmask.pdf
    for additional details on the masking. 

    return: 
        boolean array that specifies the objects that are *within* the mask 
    '''
    # same inputcount masks implemented for the imaging (this is implemented in
    # the SQL script) 
    _mask = (randoms['g_inputcount_value'] >= 4)
    _mask &= (randoms['r_inputcount_value'] >= 4)
    _mask &= (randoms['i_inputcount_value'] >= 5)
    _mask &= (randoms['z_inputcount_value'] >= 5)

    # same brightstar masks implemented for the imaging 
    _mask &= randoms['i_mask_brightstar_halo']
    _mask &= randoms['i_mask_brightstar_ghost']
    _mask &= randoms['i_mask_brightstar_blooming']

    return _mask 
