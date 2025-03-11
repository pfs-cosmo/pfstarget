--  get photometric zero-point magnitude offsets from s23b_wide.stellar_sequence_offsets
SELECT
	s1.skymap_id, 

        pq.tract, 
        pq.patch, 

        pq.ra, 
        pq.dec,

        s1.g_mag_offset, 
        s1.r_mag_offset, 
        s1.i_mag_offset, 
        s1.z_mag_offset, 
        s1.y_mag_offset
FROM
	s23b_wide.stellar_sequence_offsets as s1
        LEFT JOIN s23b_wide.patch_qa as pq USING (skymap_id) 
