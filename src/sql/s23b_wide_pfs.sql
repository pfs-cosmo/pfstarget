-- To select objects in the COSMOS field from the HSC S23B wide layer
SELECT

	-- Basic information
	f1.object_id, f1.ra, f1.dec, f1.tract, f1.patch,
	
	-- Galactic extinction for correction
	f1.a_g, f1.a_r, f1.a_i, f1.a_z, f1.a_y,
	

	-- Number of images contributing at the object center
    -- Useful to quantify the depth of the photometry
	f1.g_inputcount_value as g_input_count,
	f1.r_inputcount_value as r_input_count,
	f1.i_inputcount_value as i_input_count,
	f1.z_inputcount_value as z_input_count,
	f1.y_inputcount_value as y_input_count,

    -- Note: it seems that S23B does not isolate the photometry based on the old R2/I2 filters and provide the magnitude correction any more. I don't know whether that's because the correction has been applied or not. Might be worth checking with the HSC team.
	-- Flag for measurements from the HSC-R2/I2 filter
	-- f1.merge_measurement_r2,
	-- f1.merge_measurement_i2,
	-- c.corr_rmag,
	-- c.corr_imag,
	
	-- CModel photometry
    -- I typically download the flux and err, but it is also possible to download the magnitude and err
    -- For my case, I like to get all the CModel fluxes, including the exponential and de Vaucouluer components as a cross-check.
	-- 1. Exponential component photometry
    f1.g_cmodel_exp_flux, 
	f1.r_cmodel_exp_flux, 
	f1.i_cmodel_exp_flux, 
	f1.z_cmodel_exp_flux, 
	f1.y_cmodel_exp_flux, 
	f1.g_cmodel_exp_fluxerr as g_cmodel_exp_flux_err,
	f1.r_cmodel_exp_fluxerr as r_cmodel_exp_flux_err,
	f1.i_cmodel_exp_fluxerr as i_cmodel_exp_flux_err,
	f1.z_cmodel_exp_fluxerr as z_cmodel_exp_flux_err,
	f1.y_cmodel_exp_fluxerr as y_cmodel_exp_flux_err,
 
	-- 2. de Vaucouluer component photometry
	f1.g_cmodel_dev_flux, 
	f1.r_cmodel_dev_flux, 
	f1.i_cmodel_dev_flux, 
	f1.z_cmodel_dev_flux, 
	f1.y_cmodel_dev_flux, 
	f1.g_cmodel_dev_fluxerr as g_cmodel_dev_flux_err,
	f1.r_cmodel_dev_fluxerr as r_cmodel_dev_flux_err,
	f1.i_cmodel_dev_fluxerr as i_cmodel_dev_flux_err,
	f1.z_cmodel_dev_fluxerr as z_cmodel_dev_flux_err,
	f1.y_cmodel_dev_fluxerr as y_cmodel_dev_flux_err,
	
	-- 3. CModel photometry
	f1.g_cmodel_flux, 
	f1.r_cmodel_flux, 
	f1.i_cmodel_flux, 
	f1.z_cmodel_flux, 
	f1.y_cmodel_flux, 
	f1.g_cmodel_fluxerr as g_cmodel_flux_err,
	f1.r_cmodel_fluxerr as r_cmodel_flux_err,
	f1.i_cmodel_fluxerr as i_cmodel_flux_err,
	f1.z_cmodel_fluxerr as z_cmodel_flux_err,
	f1.y_cmodel_fluxerr as y_cmodel_flux_err,

	f1.g_cmodel_mag,
	f1.r_cmodel_mag,
	f1.i_cmodel_mag,
	f1.z_cmodel_mag,
	f1.y_cmodel_mag,
	f1.g_cmodel_magerr as g_cmodel_mag_err,
	f1.r_cmodel_magerr as r_cmodel_mag_err,
	f1.i_cmodel_magerr as i_cmodel_mag_err,
	f1.z_cmodel_magerr as z_cmodel_mag_err,
	f1.y_cmodel_magerr as y_cmodel_mag_err,

	-- 4. Flags for CModel photometry
	-- If the flag is set to True, it means the photometry is not reliable (final model fit failed)
    -- i-band flag is not included because we use it as the reference band, so we only include objects with successful i-band fits.
	f1.g_cmodel_flag,
	f1.r_cmodel_flag,
	f1.z_cmodel_flag,
	f1.y_cmodel_flag,
	
	-- PSF photometry
    -- I always rename the PSF flux columns, the original names are stupid
	f2.g_psfflux_flux as g_psf_flux, 
	f2.r_psfflux_flux as r_psf_flux, 
	f2.i_psfflux_flux as i_psf_flux, 
	f2.z_psfflux_flux as z_psf_flux, 
	f2.y_psfflux_flux as y_psf_flux, 
	f2.g_psfflux_fluxerr as g_psf_flux_err,
	f2.r_psfflux_fluxerr as r_psf_flux_err,
	f2.i_psfflux_fluxerr as i_psf_flux_err,
	f2.z_psfflux_fluxerr as z_psf_flux_err,
	f2.y_psfflux_fluxerr as y_psf_flux_err,
 
	-- PSF photometry flag (for quality cut later)
	f2.g_psfflux_flag as g_psf_flag,
	f2.r_psfflux_flag as r_psf_flag,
	f2.i_psfflux_flag as i_psf_flag,
	f2.z_psfflux_flag as z_psf_flag,
	f2.y_psfflux_flag as y_psf_flag,
	
	-- PSF-corrected aperture photometry
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diameter aperture
    -- This is the "safe choice". It is possible that a different version works the best for PFS target selection. But using a better common seeing might not work for a small fraction of the HSC footprint.
	f4.g_convolvedflux_2_15_flux as g_convolved_aper_2_15_flux, 
	f4.r_convolvedflux_2_15_flux as r_convolved_aper_2_15_flux, 
	f4.i_convolvedflux_2_15_flux as i_convolved_aper_2_15_flux, 
	f4.z_convolvedflux_2_15_flux as z_convolved_aper_2_15_flux, 
	f4.y_convolvedflux_2_15_flux as y_convolved_aper_2_15_flux, 
 
	f4.g_convolvedflux_2_15_fluxerr as g_convolved_aper_2_15_flux_err,
	f4.r_convolvedflux_2_15_fluxerr as r_convolved_aper_2_15_flux_err,
	f4.i_convolvedflux_2_15_fluxerr as i_convolved_aper_2_15_flux_err,
	f4.z_convolvedflux_2_15_fluxerr as z_convolved_aper_2_15_flux_err,
	f4.y_convolvedflux_2_15_fluxerr as y_convolved_aper_2_15_flux_err,
 
	f4.g_convolvedflux_2_15_flag as g_convolved_aper_2_15_flag,
	f4.r_convolvedflux_2_15_flag as r_convolved_aper_2_15_flag,
	f4.i_convolvedflux_2_15_flag as i_convolved_aper_2_15_flag,
	f4.z_convolvedflux_2_15_flag as z_convolved_aper_2_15_flag,
	f4.y_convolvedflux_2_15_flag as y_convolved_aper_2_15_flag,
	
	-- PSF-corrected aperture photometry **before deblending**
    -- The comparison between the deblended and undeblended fluxes can be used to identify blended objects
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diamter aperture
	f5.g_undeblended_convolvedflux_2_15_flux as g_undeblended_convolved_aper_2_15_flux, 
	f5.r_undeblended_convolvedflux_2_15_flux as r_undeblended_convolved_aper_2_15_flux, 
	f5.i_undeblended_convolvedflux_2_15_flux as i_undeblended_convolved_aper_2_15_flux, 
	f5.z_undeblended_convolvedflux_2_15_flux as z_undeblended_convolved_aper_2_15_flux, 
	f5.y_undeblended_convolvedflux_2_15_flux as y_undeblended_convolved_aper_2_15_flux, 
 
	f5.g_undeblended_convolvedflux_2_15_fluxerr as g_undeblended_convolved_aper_2_15_flux_err,
	f5.r_undeblended_convolvedflux_2_15_fluxerr as r_undeblended_convolved_aper_2_15_flux_err,
	f5.i_undeblended_convolvedflux_2_15_fluxerr as i_undeblended_convolved_aper_2_15_flux_err,
	f5.z_undeblended_convolvedflux_2_15_fluxerr as z_undeblended_convolved_aper_2_15_flux_err,
	f5.y_undeblended_convolvedflux_2_15_fluxerr as y_undeblended_convolved_aper_2_15_flux_err,
 
	f5.g_undeblended_convolvedflux_2_15_flag as g_undeblended_convolved_aper_2_15_flag,
	f5.r_undeblended_convolvedflux_2_15_flag as r_undeblended_convolved_aper_2_15_flag,
	f5.i_undeblended_convolvedflux_2_15_flag as i_undeblended_convolved_aper_2_15_flag,
	f5.z_undeblended_convolvedflux_2_15_flag as z_undeblended_convolved_aper_2_15_flag,
	f5.y_undeblended_convolvedflux_2_15_flag as y_undeblended_convolved_aper_2_15_flag,

    -- GaaP photometry: also very good for photo-z and other color-focused applications. 
    -- The _1_15x_1_5 means: Flux with 1.0 aperture after multiplying the seeing by 1.15
    f6.g_gaapflux_1_15x_1_5_flux as g_gaap_1_15x_1_5_flux,
    f6.r_gaapflux_1_15x_1_5_flux as r_gaap_1_15x_1_5_flux,
    f6.i_gaapflux_1_15x_1_5_flux as i_gaap_1_15x_1_5_flux,
    f6.z_gaapflux_1_15x_1_5_flux as z_gaap_1_15x_1_5_flux,
    f6.y_gaapflux_1_15x_1_5_flux as y_gaap_1_15x_1_5_flux,

    f6.g_gaapflux_1_15x_1_5_fluxerr as g_gaap_1_15x_1_5_flux_err,
    f6.r_gaapflux_1_15x_1_5_fluxerr as r_gaap_1_15x_1_5_flux_err,
    f6.i_gaapflux_1_15x_1_5_fluxerr as i_gaap_1_15x_1_5_flux_err,
    f6.z_gaapflux_1_15x_1_5_fluxerr as z_gaap_1_15x_1_5_flux_err,
    f6.y_gaapflux_1_15x_1_5_fluxerr as y_gaap_1_15x_1_5_flux_err,

    -- GaaP also provides an optimal aperture
    f6.g_gaapflux_1_15x_optimal_flux as g_gaap_1_15x_optimal_flux,
    f6.r_gaapflux_1_15x_optimal_flux as r_gaap_1_15x_optimal_flux,
    f6.i_gaapflux_1_15x_optimal_flux as i_gaap_1_15x_optimal_flux,
    f6.z_gaapflux_1_15x_optimal_flux as z_gaap_1_15x_optimal_flux,
    f6.y_gaapflux_1_15x_optimal_flux as y_gaap_1_15x_optimal_flux,

    f6.g_gaapflux_1_15x_optimal_fluxerr as g_gaap_1_15x_optimal_flux_err,
    f6.r_gaapflux_1_15x_optimal_fluxerr as r_gaap_1_15x_optimal_flux_err,
    f6.i_gaapflux_1_15x_optimal_fluxerr as i_gaap_1_15x_optimal_flux_err,
    f6.z_gaapflux_1_15x_optimal_fluxerr as z_gaap_1_15x_optimal_flux_err,
    f6.y_gaapflux_1_15x_optimal_fluxerr as y_gaap_1_15x_optimal_flux_err,

    -- There is also a version of GaaP run on the undeblended footprints
    f6.g_undeblended_gaapflux_1_15x_1_5_flux as g_undeblended_gaap_1_15x_1_5_flux,
    f6.r_undeblended_gaapflux_1_15x_1_5_flux as r_undeblended_gaap_1_15x_1_5_flux,
    f6.i_undeblended_gaapflux_1_15x_1_5_flux as i_undeblended_gaap_1_15x_1_5_flux,
    f6.z_undeblended_gaapflux_1_15x_1_5_flux as z_undeblended_gaap_1_15x_1_5_flux,
    f6.y_undeblended_gaapflux_1_15x_1_5_flux as y_undeblended_gaap_1_15x_1_5_flux,

    f6.g_undeblended_gaapflux_1_15x_1_5_fluxerr as g_undeblended_gaap_1_15x_1_5_flux_err,
    f6.r_undeblended_gaapflux_1_15x_1_5_fluxerr as r_undeblended_gaap_1_15x_1_5_flux_err,
    f6.i_undeblended_gaapflux_1_15x_1_5_fluxerr as i_undeblended_gaap_1_15x_1_5_flux_err,
    f6.z_undeblended_gaapflux_1_15x_1_5_fluxerr as z_undeblended_gaap_1_15x_1_5_flux_err,
    f6.y_undeblended_gaapflux_1_15x_1_5_fluxerr as y_undeblended_gaap_1_15x_1_5_flux_err,

    -- Optimal aperture for the undeblended GaaP
    f6.g_undeblended_gaapflux_1_15x_optimal_flux as g_undeblended_gaap_1_15x_optimal_flux,
    f6.r_undeblended_gaapflux_1_15x_optimal_flux as r_undeblended_gaap_1_15x_optimal_flux,
    f6.i_undeblended_gaapflux_1_15x_optimal_flux as i_undeblended_gaap_1_15x_optimal_flux,
    f6.z_undeblended_gaapflux_1_15x_optimal_flux as z_undeblended_gaap_1_15x_optimal_flux,
    f6.y_undeblended_gaapflux_1_15x_optimal_flux as y_undeblended_gaap_1_15x_optimal_flux,

    f6.g_undeblended_gaapflux_1_15x_optimal_fluxerr as g_undeblended_gaap_1_15x_optimal_flux_err,
    f6.r_undeblended_gaapflux_1_15x_optimal_fluxerr as r_undeblended_gaap_1_15x_optimal_flux_err,
    f6.i_undeblended_gaapflux_1_15x_optimal_fluxerr as i_undeblended_gaap_1_15x_optimal_flux_err,
    f6.z_undeblended_gaapflux_1_15x_optimal_fluxerr as z_undeblended_gaap_1_15x_optimal_flux_err,
    f6.y_undeblended_gaapflux_1_15x_optimal_fluxerr as y_undeblended_gaap_1_15x_optimal_flux_err,
	
	-- SDSS Shape without PSF correction (Using i-band; can use others too)
	f2.i_sdssshape_shape11 as i_sdss_shape_11,
	f2.i_sdssshape_shape12 as i_sdss_shape_12,
	f2.i_sdssshape_shape22 as i_sdss_shape_22,
	f2.i_sdssshape_shape11err as i_sdss_shape_11_err,
	f2.i_sdssshape_shape12err as i_sdss_shape_12_err,
	f2.i_sdssshape_shape22err as i_sdss_shape_22_err,
	
	-- Shape of the CModel model (i & r-band)
    -- This takes the PSF convolution into account and one can recover the size, ellipticity, and position angle of the object from these three moments. 
    -- I like to use the exponential fit one because it is more stable.
	m1.i_cmodel_exp_ellipse_11, 
	m1.i_cmodel_exp_ellipse_22, 
	m1.i_cmodel_exp_ellipse_12,
	m1.i_cmodel_ellipse_11, 
	m1.i_cmodel_ellipse_22, 
	m1.i_cmodel_ellipse_12,
	
	m1.r_cmodel_exp_ellipse_11, 
	m1.r_cmodel_exp_ellipse_22, 
	m1.r_cmodel_exp_ellipse_12,
	m1.r_cmodel_ellipse_11, 
	m1.r_cmodel_ellipse_22, 
	m1.r_cmodel_ellipse_12,

	-- Extendedness of the object
	f1.g_extendedness_value,
	f1.r_extendedness_value,
	f1.i_extendedness_value,
	f1.z_extendedness_value,
	f1.y_extendedness_value,
	f1.g_extendedness_flag,
	f1.r_extendedness_flag,
	f1.i_extendedness_flag,
	f1.z_extendedness_flag,
	f1.y_extendedness_flag,

    -- Blendedness and deblender-related diagnoses.
    m2.g_deblend_blendedness as g_blendedness,
    m2.r_deblend_blendedness as r_blendedness,
    m2.i_deblend_blendedness as i_blendedness,
    m2.z_deblend_blendedness as z_blendedness,
    m2.y_deblend_blendedness as y_blendedness,

    -- This is the "de-noised" version of the blendedness (see Bosch et al. 2018)
    m2.g_blendedness_abs,
    m2.r_blendedness_abs,
    m2.i_blendedness_abs,
    m2.z_blendedness_abs,
    m2.y_blendedness_abs,

    -- The maxoverlap metric: Maximum overlap with all of the other neighbors flux combined (not familiar with this, but sounds useful).
    m2.g_deblend_maxoverlap, 
    m2.r_deblend_maxoverlap,
    m2.i_deblend_maxoverlap,
    m2.z_deblend_maxoverlap,
    m2.y_deblend_maxoverlap,

    -- The deblender flag
    m2.g_blendedness_flag,
    m2.r_blendedness_flag,
    m2.i_blendedness_flag,
    m2.z_blendedness_flag,
    m2.y_blendedness_flag,
	
	-- Flags for later selection
	-- 1. The general failure flag
	f1.g_pixelflags,
	f1.r_pixelflags,
	f1.i_pixelflags,
	f1.z_pixelflags,
	f1.y_pixelflags,
	
	-- 2. Saturated or interpolated pixels on the footprint (not center)
	f1.g_pixelflags_saturated,
	f1.r_pixelflags_saturated,
	f1.i_pixelflags_saturated,
	f1.z_pixelflags_saturated,
	f1.y_pixelflags_saturated,
	f1.g_pixelflags_interpolated,
	f1.r_pixelflags_interpolated,
	f1.i_pixelflags_interpolated,
	f1.z_pixelflags_interpolated,
	f1.y_pixelflags_interpolated,

	-- 3. Other pixel flags
	f1.g_pixelflags_bad,
	f1.r_pixelflags_bad,
	f1.i_pixelflags_bad,
	f1.z_pixelflags_bad,
	f1.y_pixelflags_bad,
	f1.g_pixelflags_suspectcenter,
	f1.r_pixelflags_suspectcenter,
	f1.i_pixelflags_suspectcenter,
	f1.z_pixelflags_suspectcenter,
	f1.y_pixelflags_suspectcenter,
	f1.g_pixelflags_clippedcenter,
	f1.r_pixelflags_clippedcenter,
	f1.i_pixelflags_clippedcenter,
	f1.z_pixelflags_clippedcenter,
	f1.y_pixelflags_clippedcenter,

	-- 4. Bright object masks
	f1.g_pixelflags_bright_object,
	f1.r_pixelflags_bright_object,
	f1.i_pixelflags_bright_object,
	f1.z_pixelflags_bright_object,
	f1.y_pixelflags_bright_object,

    -- 5. Detailed bright star masks 
    -- It is also possible to gather more detailed information about the bright star masks here.
    msk.g_mask_brightstar_halo, 
    msk.g_mask_brightstar_dip, 
    msk.g_mask_brightstar_ghost,
    msk.g_mask_brightstar_blooming,
    msk.g_mask_brightstar_ghost12,
    msk.g_mask_brightstar_ghost15,

    msk.r_mask_brightstar_halo,
    msk.r_mask_brightstar_dip,
    msk.r_mask_brightstar_ghost,
    msk.r_mask_brightstar_blooming,
    msk.r_mask_brightstar_ghost12,
    msk.r_mask_brightstar_ghost15,

    msk.i_mask_brightstar_halo,
    msk.i_mask_brightstar_dip,
    msk.i_mask_brightstar_ghost,
    msk.i_mask_brightstar_blooming,
    msk.i_mask_brightstar_ghost12,
    msk.i_mask_brightstar_ghost15,

    msk.z_mask_brightstar_halo,
    msk.z_mask_brightstar_dip,
    msk.z_mask_brightstar_ghost,
    msk.z_mask_brightstar_blooming,
    msk.z_mask_brightstar_ghost12,
    msk.z_mask_brightstar_ghost15

	msk.y_mask_brightstar_halo,
    msk.y_mask_brightstar_dip,
    msk.y_mask_brightstar_ghost,
    msk.y_mask_brightstar_blooming,
    msk.y_mask_brightstar_ghost12,
    msk.y_mask_brightstar_ghost15
	
    -- I don't think the S23B photo-z is ready yet, so I will skip this part for now.
	-- Mizuki photo-z information 
	-- p1.photoz_mean as pz_mean_mizuki, 
	-- p1.photoz_best as pz_best_mizuki, 
	-- p1.photoz_conf_mean as pz_conf_mean_mizuki, 
	-- p1.photoz_conf_best as pz_conf_best_mizuki, 
	-- p1.photoz_risk_mean as pz_risk_mean_mizuki, 
	-- p1.photoz_risk_best as pz_risk_best_mizuki,
	-- p1.photoz_std_mean as pz_std_mean_mizuki, 
	-- p1.photoz_std_best as pz_std_best_mizuki,
	-- p1.photoz_err68_min as pz_err68_min_mizuki, 
	-- p1.photoz_err68_max as pz_err68_max_mizuki, 
	-- p1.stellar_mass as mstar_mizuki,
	-- p1.stellar_mass_err68_min as mstar_min_mizuki,
	-- p1.stellar_mass_err68_max as mstar_max_mizuki,
	
	-- DNNZ photo-z information 
	-- p2.photoz_mean as pz_mean_dnnz, 
	-- p2.photoz_best as pz_best_dnnz, 
	-- p2.photoz_conf_mean as pz_conf_mean_dnnz, 
	-- p2.photoz_conf_best as pz_conf_best_dnnz, 
	-- p2.photoz_risk_mean as pz_risk_mean_dnnz, 
	-- p2.photoz_risk_best as pz_risk_best_dnnz,
	-- p2.photoz_std_mean as pz_std_mean_dnnz, 
	-- p2.photoz_std_best as pz_std_best_dnnz,
	-- p2.photoz_err68_min as pz_err68_min_dnnz, 
	-- p2.photoz_err68_max as pz_err68_max_dnnz, 	

	-- DEMP photo-z information 
	-- p3.photoz_mean as pz_mean_demp, 
	-- p3.photoz_best as pz_best_demp, 
	-- p3.photoz_conf_mean as pz_conf_mean_demp, 
	-- p3.photoz_conf_best as pz_conf_best_demp, 
	-- p3.photoz_risk_mean as pz_risk_mean_demp, 
	-- p3.photoz_risk_best as pz_risk_best_demp,
	-- p3.photoz_std_mean as pz_std_mean_demp, 
	-- p3.photoz_std_best as pz_std_best_demp,
	-- p3.photoz_err68_min as pz_err68_min_demp, 
	-- p3.photoz_err68_max as pz_err68_max_demp	

FROM
	s23b_wide.forced as f1
	LEFT JOIN s23b_wide.forced2 as f2 USING (object_id)
	LEFT JOIN s23b_wide.forced4 as f4 USING (object_id)
	LEFT JOIN s23b_wide.forced5 as f5 USING (object_id)
	LEFT JOIN s23b_wide.forced6 as f6 USING (object_id)
	LEFT JOIN s23b_wide.meas as m1 USING (object_id)
	LEFT JOIN s23b_wide.meas2 as m2 USING (object_id)
    LEFT JOIN s23b_wide.masks as msk USING (object_id)
    -- Mag correction is not available in S23B
	-- LEFT JOIN s23b_wide.mag_corr as c USING (object_id)
    -- Photo-z not available in S23B
	-- LEFT JOIN s23b_wide.photoz_mizuki as p1 USING (object_id)
	-- LEFT JOIN s23b_wide.photoz_dnnz as p2 USING (object_id)
	-- LEFT JOIN s23b_wide.photoz_demp as p3 USING (object_id)

WHERE
       f1.isprimary
	-- Region
	-- This is for the COSMOS field
	-- AND boxSearch(coord, 149.22484, 150.81303, 1.541319, 2.87744)

    -- Tract 
	AND m1.tract IN ({$tract})

	-- Photometric cut
	-- Note: will do mag cut in post-processing
	-- AND f1.i_cmodel_mag <= 23.5
	AND NOT f1.i_cmodel_flag

	-- Extendedness cut 
	-- Note: this is done in target seleciton post-processing
    -- Note: need to think about this cut. Is it true that all potential PSF targets are extended objects?
	-- AND f1.r_extendedness_value > 0
	-- AND f1.i_extendedness_value > 0

	-- Full-depth full-color cut
	AND f1.g_inputcount_value >= 4
	AND f1.r_inputcount_value >= 4
	AND f1.i_inputcount_value >= 5
	AND f1.z_inputcount_value >= 5
	AND f1.y_inputcount_value >= 5

	-- Not failed at finding the center
	AND NOT f2.g_sdsscentroid_flag
	AND NOT f2.r_sdsscentroid_flag
	AND NOT f2.i_sdsscentroid_flag
	AND NOT f2.z_sdsscentroid_flag
	AND NOT f2.y_sdsscentroid_flag

	-- The object's center is not outside the image
	AND NOT f1.g_pixelflags_edge    
	AND NOT f1.r_pixelflags_edge
	AND NOT f1.i_pixelflags_edge
	AND NOT f1.z_pixelflags_edge
	AND NOT f1.y_pixelflags_edge

	-- Not saturated at the center  
	AND NOT f1.g_pixelflags_saturatedcenter
    AND NOT f1.r_pixelflags_saturatedcenter
    AND NOT f1.i_pixelflags_saturatedcenter
    AND NOT f1.z_pixelflags_saturatedcenter
    AND NOT f1.y_pixelflags_saturatedcenter

	-- The center is not interpolated
	AND NOT f1.g_pixelflags_interpolatedcenter
	AND NOT f1.r_pixelflags_interpolatedcenter
	AND NOT f1.i_pixelflags_interpolatedcenter
	AND NOT f1.z_pixelflags_interpolatedcenter
	AND NOT f1.y_pixelflags_interpolatedcenter

	-- The center is not affected by a cosmic ray
	AND NOT f1.g_pixelflags_crcenter
    AND NOT f1.r_pixelflags_crcenter
    AND NOT f1.i_pixelflags_crcenter
    AND NOT f1.z_pixelflags_crcenter
    AND NOT f1.y_pixelflags_crcenter
;