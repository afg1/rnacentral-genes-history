#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fetch_ensembl_prefixes {
    input:
        path(taxon_query)

    output:
        path('taxid_name.csv')

    script:
    """
    psql -v ON_ERROR_STOP=1 -f ${taxon_query} "$PGDATABASE" > taxid_name.csv
    """
}

process fetch_regions_data {
    maxForks 50
    input:
        tuple val(taxid), path(regions_query)
    output:
        tuple val(taxid), path('regions.csv')
    script:
    """
    psql -v ON_ERROR_STOP=1 -v taxid=${taxid} -f ${regions_query} "$PGDATABASE" > regions.csv
    """
}

process get_organism_paths {

    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    
    input:
        tuple path(raw_taxid_name), path(script_path)

    output:
        path('taxid_name_dirname.csv')
    
    script:
    """
    python ${script_path} ${raw_taxid_name} taxid_name_dirname.csv
    """
}

process fetch_inactive_ids {
    tag "${release}: inactive ids copy"
    queue 'datamover'
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    input:
        tuple val(release), val(dummy)
    
    output:
        tuple val(release), path("inactive_ids_release_${release}")
    
    script:
    """
    cp /nfs/ftp/public/databases/RNAcentral/releases/${release}.0/sequences/rnacentral_inactive.fasta.gz .
    gzip -d rnacentral_inactive.fasta.gz
    grep -oE 'URS[0-9A-Z]+' rnacentral_inactive.fasta > inactive_ids_release_${release}
    rm rnacentral_inactive.fasta
    """
}

process copy_gff {
    tag "Release ${meta.release}: ${meta.org_name} GFF copy"
    queue 'datamover'
    errorStrategy 'ignore'

    input:
        tuple val(meta), val(dummy)

    output:
        tuple val(meta), path('*.gff3')

    script:
    """
    cp /nfs/ftp/public/databases/RNAcentral/releases/${meta.release}.0/genome_coordinates/gff3/${meta.dirname}.*.gff3.gz .
    gzip -d ${meta.dirname}.*.gff3.gz
    """
}

process convert_gff_to_parquet {
    tag "${meta.release}: ${meta.org_name} GFF conv"
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    input:
        tuple val(meta), path(input_gff), path(regions_file)
    
    output:
        tuple val(meta), path('*.parquet'), path(regions_file)

    script:
    """
    rnac genes utils convert --gff_file ${input_gff} --taxid ${meta.taxid} --regions_data ${regions_file}
    """
}

process preprocess_transcripts {
    tag "Release ${meta.release}: ${meta.org_name} preprocessing"
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    memory { 256.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(meta), path(input_parquet), path(regions_file),  path(so_model)
    
    output:
        tuple val(meta), path("${meta.dirname}_features.parquet")

    script:
    """
    rnac genes infer preprocess \
    --transcripts_file ${input_parquet} \
    --so_model_path ${so_model} \
    --regions_data ${regions_file} \
    --output ${meta.dirname}_features.parquet \
    --no-parallel
    """

}

process classify_pairs {
    tag "Release ${meta.release}: ${meta.org_name} classification"
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    memory { 1.GB * task.attempt }
    errorStrategy 'ignore' //{ task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1
    cpus 4


    input:
        tuple val(meta), path(transcripts), path(regions), path(features), path(rf_model)

    output:
        tuple val(meta), path("*.json")

    script:
    """
    export OMP_NUM_THREADS=4

    rnac genes infer classify \
    --transcripts_file ${transcripts} \
    --features_file ${features} \
    --taxid ${meta.taxid} \
    --output_dir . \
    --model_path ${rf_model}

    mv genes_${meta.taxid}.json genes_${meta.taxid}_${meta.release}.json
    """
}

process forward_merge {
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    tag "Forward_merge ${taxid}"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    cpus 4

    input:
        tuple val(taxid), path(genes_files), path(inactive_files), val(releases)

    output:
        tuple val(taxid), path("final_merged_${taxid}.json")

    script:
    """
    # Create associative arrays for genes and inactive files by release
    declare -A genes_files
    declare -A inactive_files

    # Parse genes files: genes_TAXID_RELEASE.json
    for file in ${genes_files}; do
        if [[ \$file =~ genes_${taxid}_([0-9]+)\\.json ]]; then
            release=\${BASH_REMATCH[1]}
            genes_files[\$release]=\$file
        fi
    done

    # Parse inactive files: inactive_ids_release_RELEASE
    for file in ${inactive_files}; do
        if [[ \$file =~ inactive_ids_release_([0-9]+) ]]; then
            release=\${BASH_REMATCH[1]}
            inactive_files[\$release]=\$file
        fi
    done

    # Use the provided releases list (already sorted)
    releases=(${releases.join(' ')})
    echo "Processing releases in order: \${releases[*]}"

    # Find the first available release for this taxon
    first_available_release=""
    for release in "\${releases[@]}"; do
        if [[ -n "\${genes_files[\$release]}" && -n "\${inactive_files[\$release]}" ]]; then
            first_available_release=\$release
            echo "First available release for taxon ${taxid}: \$first_available_release"
            break
        fi
    done

    # Check if we found any files at all
    if [[ -z "\$first_available_release" ]]; then
        echo "ERROR: No genes or inactive files found for taxon ${taxid} in any release"
        exit 1
    fi

    # Start with the first available release
    prev_file="\${genes_files[\$first_available_release]}"
    prev_release="\$first_available_release"

    echo "Starting with release \$first_available_release"

    # Loop through remaining releases and merge sequentially  
    for release in "\${releases[@]}"; do
        # Skip if this release is before or equal to our starting release
        if [[ \$release -le \$first_available_release ]]; then
            continue
        fi
        
        curr_file="\${genes_files[\$release]}"
        inactive_file="\${inactive_files[\$release]}"
        
        # Skip this release if files are missing
        if [[ -z "\$curr_file" || -z "\$inactive_file" ]]; then
            echo "Skipping release \$release: missing genes_file='\$curr_file' or inactive_file='\$inactive_file'"
            continue
        fi
        
        output_file="\${release}_merged.json"
        
        echo "Merging release \$prev_release + \$release"
        echo "  Previous: \$prev_file"
        echo "  Current: \$curr_file"
        echo "  Inactive: \$inactive_file"
        echo "  Output: \$output_file"
        
        rnac genes merge \\
            --previous_genes "\$prev_file" \\
            --next_genes "\$curr_file" \\
            --output "\$output_file" \\
            --inactive_ids "\$inactive_file" \\
            --prev_release_number \$prev_release \\
            --next_release_number \$release
        
        # Update for next iteration
        prev_file="\$output_file"
        prev_release=\$release
    done

    # Copy final result to expected output name
    cp "\$prev_file" "final_merged_${taxid}.json"
    """
}


process calculate_metadata {
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    tag "${taxid} metadata calculation"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    cpus 4
    maxForks 10

    input:
        tuple val(taxid), path(merged_genes)

    output:
        tuple val(taxid), path("${taxid}_metadata.json")

    script:
    """
    rnac genes utils process-metadata ${merged_genes} ${taxid}_metadata.json
    """
}

process upload_genes {
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    tag "${taxid} gene upload"
    maxForks 1

    input:
        tuple val(taxid), path(merged_genes)

    output:
        tuple val(taxid), val("DONE GENES")

    script:
    """
    rnac genes utils store-genes ${merged_genes} --taxid ${taxid}
    """
}

process upload_metadata {
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    tag "${taxid} metadata upload"
    maxForks 1

    input:
        tuple val(taxid), path(metadata)

    output:
        tuple val(taxid), val("DONE METADATA")

    script:
    """
    rnac genes utils store-metadata ${metadata}
    """
}


workflow {
    taxa_query = Channel.fromPath('/hps/nobackup/agb/rnacentral/rnacentral-genes-history/sql/get_taxids.sql')
    release_file = Channel.fromPath('/hps/nobackup/agb/rnacentral/rnacentral-genes-history/releases.txt')
    so_model = Channel.fromPath('/hps/nobackup/agb/rnacentral/rnacentral-genes-history/so_model.emb')
    rf_model = Channel.fromPath('/hps/nobackup/agb/rnacentral/rnacentral-genes-history/rf_model.onnx')
    org_name_script = Channel.fromPath("/hps/nobackup/agb/rnacentral/rnacentral-genes-history/utils/organism_name.py")
    regions_query = Channel.fromPath('/hps/nobackup/agb/rnacentral/rnacentral-genes-history/sql/dump_regions.sql')

    taxid_name_dirname = taxa_query | fetch_ensembl_prefixes | combine(org_name_script) | get_organism_paths | splitCsv(header: true) | map { row -> [row.taxid, row.organism_name, row.transformed_name] }
    releases = release_file | splitText | map { it.trim() } | filter { it != "" }  
    releases_list = releases.toSortedList { a, b -> a as Integer <=> b as Integer }
    unique_taxids = taxid_name_dirname
        .map { taxid, org_name, dirname -> taxid }
        .unique()
    
    // Fetch regions data for each taxid
    regions_data = unique_taxids
        .combine(regions_query)
        .map { taxid, query_file -> [taxid, query_file] }
        | fetch_regions_data


    combo = taxid_name_dirname.combine(releases).map { taxid, org_name, dirname, release ->
            // Create a meta map containing ALL the metadata we want to track
            def meta = [
                taxid: taxid,
                org_name: org_name, 
                dirname: dirname,
                release: release,
                sample_id: "${org_name}_rel${release}"  // Computed field
            ]
            
            // Return tuple of [meta, files]
            // The empty list [] is for initial files (none in this case)
            [meta, []]
        }

    meta_by_release = combo
        .map { meta, dummy -> [meta.release, meta] }
        .groupTuple()  // Groups all taxa for each release

    inactive_ids = releases.map { release -> [release, []] }  // Convert to [release, dummy] structure
    | fetch_inactive_ids

    gff_files = combo | copy_gff
    gff_with_regions = gff_files
        .map { meta, gff_file -> [meta.taxid, meta, gff_file] }
        .combine(regions_data, by: 0)
        .map { taxid, meta, gff_file, regions_file -> [meta, gff_file, regions_file] }

    transcripts_with_regions = gff_with_regions | convert_gff_to_parquet

    features = preprocess_transcripts(transcripts_with_regions.combine(so_model))

    genes = classify_pairs(transcripts_with_regions.join(features).combine(rf_model))
    
    genes_collected = genes
        .map { meta, json_file -> 
            [meta.taxid, json_file]
        }
        .groupTuple()
    
    inactive_with_meta = inactive_ids
    .join(meta_by_release)  // Join on release number
    .flatMap { release, inactive_file, meta_list ->
        // Create one entry per taxon for this release
        meta_list.collect { meta -> [meta, inactive_file] }
    }


    inactive_collected = inactive_with_meta
        .map { meta, inactive_file -> [meta.taxid, inactive_file] }
        .groupTuple()

    combined_for_merge = genes_collected.join(inactive_collected).combine(releases_list)

    merged_genes = combined_for_merge | forward_merge

    gene_metadata = merged_genes | calculate_metadata

    // genes_done = merged_genes | upload_genes

    // metadata_done = gene_metadata | upload_metadata

}