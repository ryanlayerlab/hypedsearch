
/*SELECT * 
FROM kmers 
where mass between 226 and 227
order by protein, location_start*/

--SELECT *
--FROM proteins
--where id = 30


/*SELECT *
FROM kmers
where mass = 261.6369194267578*/

Select *
FROM kmers
where protein=30
and location_start = 13
and location_end = 16
and ion = 1
and charge = 2
