-- To select objects from the HSC S23B wide layer
SELECT *
FROM
	s23b_wide.random as r1 
	LEFT JOIN s23b_wide.random_masks as r2 USING (object_id)
WHERE
    -- Tract 
    r1.tract IN ({$tract})
