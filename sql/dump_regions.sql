COPY(
  SELECT 
   json_build_object(
        'region_name', sr.region_name,
        'assembly_id', sr.assembly_id,
        'region_id', sr.id,
    )
  sr.assembly_id as assembly_id,
  sr.id as region_id

  FROM rnc_sequence_regions sr 

  WHERE sr.urs_taxid LIKE '%_' || :'taxid'
) TO STDOUT;
