'''

module for correcting for galactic extinction


'''

def _extinction_correct(hsc, method='sfd98'): 
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
        
        # get photometry zero-point offsets from wide.stellar_sequence_offset
        # (see Issue #8 for details)  
        grizy_offset = _get_zeropoint_correct(hsc['tract'], hsc['patch']) 

        # correct magnitudes for dust extinction
        g_mag = hsc["g_cmodel_mag"] - g_a - grizy_offset[0]
        r_mag = hsc["r_cmodel_mag"] - r_a - grizy_offset[1]
        i_mag = hsc["i_cmodel_mag"] - i_a - grizy_offset[2]
        z_mag = hsc["z_cmodel_mag"] - z_a - grizy_offset[3]
        y_mag = hsc["y_cmodel_mag"] - y_a - grizy_offset[4]

    else:
        # IMPLEMENT RONGPU'S DUST EXTINCTION CORRECTION HERE 
        raise NotImplementedError

    return g_mag, r_mag, i_mag, z_mag, y_mag 


def _get_zeropoint_correct(tract, patch): 
    ''' return g/r/i/z/y-band photometric zeropoint correction based on tract
    and patch. 
    '''
    # implement look up table for pdr3_wide.stellar_sequence_offsets based on
    # tract and patch

    return offset 

