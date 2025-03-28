-- To select objects from the HSC S23B wide layer
SELECT *
FROM
	s23b_wide.random as r1 
	LEFT JOIN s23b_wide.random_masks as r2 USING (object_id)
WHERE
    -- Tract 
    r1.tract IN ({$tract})
    
    ----- same inputcount masks implemented for the imaging 
    ---AND r1.g_inputcount_value >= 4
    ---AND r1.r_inputcount_value >= 4
    ---AND r1.i_inputcount_value >= 5
    ---AND r1.z_inputcount_value >= 5
