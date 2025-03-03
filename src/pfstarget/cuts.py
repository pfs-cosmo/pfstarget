'''

module for target selection 


'''
import numpy as np 



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

    # quality cuts
    is_quality = quality_cuts(objects)

    # star-galaxy separation
    is_galaxy = star_galaxy(objects, cut=star_galaxy_cut)

    # color cut
    is_color = color_cut(objects, magnitude_cut=magnitude_cut) 

    return is_quality & is_galaxy & is_color 


def color_cut(objects, magnitude_cut=22.5, g_r_cut=0.15, color_slope=2.0, color_yint=-0.15): 
    ''' impose color cut to select ELG within the redshift range of 0.6 < z <
    2.4
    '''
    # faint magnitude cut 
    cuts = objects['G_MAG'] > magnitude_cut  
    
    # color cut 
    cuts &= (((objects['G_MAG'] - objects['R_MAG']) < 0.15) | # g-r cut (for 1.6 < z < 2.4 ELGs) 
             ((objects['I_MAG'] - objects['Y_MAG']) 
              > color_slope * (objects['G_MAG'] - objects['R_MAG'])- color_yint))
    return cuts 


def quality_cuts(objects): 
    ''' quality cuts. Currently we only impose a magnitude depend cut on g-band
    magnitude error.
    '''
    # g-band magnitude error cut
    cuts = (objects['G_ERR'] < objects['G_MAG'] * 0.05 - 1.1) 

    # TODO: extreme colors
    # There should probably be a cut on extreme colors 
    return cuts


def star_galaxy(i_cmodel, i_psf, cut=-0.15): 
    ''' star-galaxy separation. Currently PSF star-galaxy separation is:

    [i-band cmodel magnitude] - [i-band PSF magnitude] < cut

    This cut essentially determines whether the object is extended and has
    significant light outside of the PSF profile. 
    '''
    return (i_cmodel - i_psf < cut)


def _prepare_hsc(hsc, dust_extinction='default'): 
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
             ('I_MEAS_PSF_MAG', 'f4')
             ]

    objects = np.zeros(len(hsc), dtype=dtype)

    # grizy magnitudes corrections for galactic dust extinction 
    if dust_extinction == 'default': 
        # E(B -V)
        # https://hsc-release.mtk.nao.ac.jp/schema/#pdr3.pdr3_wide.forced
        # 3.24::real * "_forced:part1".extinction_bv AS a_g,
        # 2.276::real * "_forced:part1".extinction_bv AS a_r,
        # 1.633::real * "_forced:part1".extinction_bv AS a_i,
        # 1.263::real * "_forced:part1".extinction_bv AS a_z,
        # 1.075::real * "_forced:part1".extinction_bv AS a_y,
        g_a = hsc["a_g"]
        r_a = hsc["a_r"]
        i_a = hsc["a_i"]
        z_a = hsc["a_z"]
        y_a = hsc["a_y"]
    else: 
        # IMPLEMENT RONGPU'S DUST EXTINCTION CORRECTION HERE 
        raise NotImplementedError

    # correct magnitudes for dust extinction
    objects['G_MAG'] = hsc["g_cmodel_mag"] - g_a
    objects['R_MAG'] = hsc["r_cmodel_mag"] - r_a
    objects['I_MAG'] = hsc["i_cmodel_mag"] - i_a
    objects['Z_MAG'] = hsc["z_cmodel_mag"] - z_a
    objects['Y_MAG'] = hsc["y_cmodel_mag"] - y_a
           
    # g-band CMODEL magnitude used for quality cuts 
    objects['G_ERR'] = hsc["g_cmodel_mag_err"]
    objects['R_ERR'] = hsc["r_cmodel_mag_err"]
    objects['I_ERR'] = hsc["i_cmodel_mag_err"]
    objects['Z_ERR'] = hsc["z_cmodel_mag_err"]
    objects['Y_ERR'] = hsc["y_cmodel_mag_err"]

    # i-band measured CMODEL and PSF magnitudes for star-galaxy separation 
    objects['I_PSF_MAG'] = hsc["i_psf_mag"]
    # commented out the lines below because the updated HSC tract data no
    # longer includes measured cmodel or PSF magnitudes (see issue #4). 
    #objects['I_MEAS_CMODEL_MAG'] = hsc["meas_i_cmodel_mag"]
    #objects['I_MEAS_PSF_MAG'] = hsc["meas_i_psfflux_mag"]

    # additional columns 
    objects['OBJID'] = hsc['object_id']  # object id 
    objects['RA'] = hsc['i_ra'] # RA from i-band measurement
    objects['DEC'] = hsc['i_dec'] # RA from i-band measurement

    return objects
