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
    tag "${meta.release}: ${meta.org_name} inactive ids copy"
    queue 'datamover'

    input:
        tuple val(meta), val(dummy)
    
    output:
        tuple val(meta), path("inactive_ids_release_${meta.release}")
    
    script:
    """
    cp /nfs/ftp/public/databases/RNAcentral/releases/${meta.release}.0/sequences/rnacentral_inactive.fasta.gz .
    gzip -d rnacentral_inactive.fasta.gz
    grep -oE 'URS[0-9A-Z]+' rnacentral_inactive.fasta > inactive_ids_release_${meta.release}
    rm rnacentral_inactive.fasta
    """
}

process copy_gff {
    tag "${meta.release}: ${meta.org_name} GFF copy"
    queue 'datamover'

    input:
        tuple val(meta), val(dummy)

    output:
        tuple val(meta), path('*.gff')

    script:
    """
    cp /nfs/ftp/public/databases/RNAcentral/releases/${meta.release}.0/genome_coordinates/gff3/${meta.dirname}.gff3.gz .
    gzip -d ${meta.dirname}.gff3.gz
    """
}

process convert_gff_to_parquet {
    tag "${meta.release}: ${meta.org_name} GFF conv"
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'

    input:
        tuple val(meta), path(input_gff)
    
    output:
        tuple val(meta), path('*.parquet')

    script:
    """
    rnac genes utils convert --gff_file ${input_gff} --taxid ${meta.taxid}
    """
}

process preprocess_transcripts {
    tag "${meta.release}: ${meta.org_name} preprocessing"
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'

    input:
        tuple val(meta), path(input_parquet), path(so_model)
    
    output:
        tuple val(meta), path("${meta.dirname}_features.parquet")

    script:
    """
    rnac genes infer preprocess \
    --transcripts_file ${input_parquet} \
    --so_model_path ${so_model} \
    --output ${meta.dirname}_features.parquet \
    --no-parallel
    """

}

process classify_pairs {
    tag "${meta.release}: ${meta.org_name} classification"

    input:
        tuple val(meta), path(transcripts), path(features), path(rf_model)

    output:
        tuple val(meta), path("*.json")

    script:
    """
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
    
    # Debug: show what files we found
    for release in "\${releases[@]}"; do
        echo "Release \$release: genes=\${genes_files[\$release]}, inactive=\${inactive_files[\$release]}"
    done
    
    # Start with the first release
    prev_file="\${genes_files[\${releases[0]}]}"
    prev_release="\${releases[0]}"
    
    # Loop through remaining releases and merge sequentially  
    for i in "\${!releases[@]}"; do
        if [ \$i -eq 0 ]; then
            continue  # Skip first release
        fi
        
        curr_release="\${releases[\$i]}"
        curr_file="\${genes_files[\$curr_release]}"
        inactive_file="\${inactive_files[\$curr_release]}"
        output_file="\${curr_release}_merged.json"
        
        echo "Merging release \$prev_release + \$curr_release"
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
            --next_release_number \$curr_release
        
        # Update for next iteration
        prev_file="\$output_file"
        prev_release=\$curr_release
    done
    
    # Copy final result to expected output name
    cp "\$prev_file" "final_merged_${taxid}.json"
    """
}


process calculate_metadata {
    container 'oras://ghcr.io/rnacentral/rnacentral-import-pipeline:latest'
    tag "${taxid} metadata calculation"

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
    tag "${taxid} metadata calculation"
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
    tag "${taxid} metadata calculation"
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
    taxa_query = Channel.fromPath('sql/get_taxids.sql')
    release_file = Channel.fromPath('releases.txt')
    so_model = Channel.fromPath('so_model.emb')
    rf_model = Channel.fromPath('rf_model.onnx')
    org_name_script = Channel.fromPath("utils/organism_name.py")

    taxid_name_dirname = taxa_query | fetch_ensembl_prefixes | combine(org_name_script) | get_organism_paths | splitCsv(header: true) | map { row -> [row.taxid, row.org_name, row.transformed_name] }
    releases = release_file | splitText | map { it.trim() } | filter { it != "" }  
    releases_list = releases.collect().map { it.sort { a, b -> a as Integer <=> b as Integer } }

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

    inactive_ids = combo | fetch_inactive_ids

    transcripts = combo | copy_gff | convert_gff_to_parquet 

    features = preprocess_transcripts(transcripts.combine(so_model) )

    genes = classify_pairs(transcripts.join(features).combine(rf_model))
    
    genes_collected = genes
        .map { meta, json_file -> 
            [taxid, json_file]
        }
        .groupTuple()
    
    inactive_collected = inactive_ids
        .map { meta, inactive_file ->
            [taxid, inactive_file]
        }
        .groupTuple()

    combined_for_merge = genes_collected.join(inactive_collected).combine(releases_list)

    merged_genes = combined_for_merge | forward_merge

    gene_metadata = merged_genes | calculate_metadata

    genes_done = merged_genes | upload_genes

    metadata_done = gene_metadata | upload_metadata

}