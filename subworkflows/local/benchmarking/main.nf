/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BENCHMARKING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Compares taxonomic profiles against a gold standard.
    Produces a JSON metrics file consumed by the stimulus optimization framework.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CONVERT_TO_BIOBOXES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(taxpasta_tsv)
    path taxdump_dir

    output:
    tuple val(meta), path("*.bioboxes.profile"), emit: profile
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert_to_bioboxes.py \\
        --input ${taxpasta_tsv} \\
        --sample-id ${meta.id} \\
        --taxdump ${taxdump_dir} \\
        --output ${prefix}.bioboxes.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        convert_to_bioboxes: 1.0.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_PROFILE > ${prefix}.bioboxes.profile
    @Version:0.9.1
    @SampleID:${meta.id}
    @Ranks:superkingdom|phylum|class|order|family|genus|species|strain
    @@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE
    END_PROFILE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        convert_to_bioboxes: 1.0.0
    END_VERSIONS
    """
}

process EVALUATE_PROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(query_profile)
    path truth_profile
    val tools_used

    output:
    tuple val(meta), path("metrics.json"), emit: metrics
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    evaluate_profile.py \\
        --prediction ${query_profile} \\
        --truth ${truth_profile} \\
        --sample ${meta.id} \\
        --tools "${tools_used}" \\
        --output metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        evaluate_profile: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_JSON > metrics.json
    {
        "sample": "${meta.id}",
        "tool_combination": "${tools_used}",
        "per_rank_metrics": {},
        "summary": {
            "weighted_f1": 0.0,
            "species_precision": 0.0,
            "species_recall": 0.0,
            "species_f1": 0.0,
            "objective_rank": "species"
        }
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        evaluate_profile: 1.0.0
    END_VERSIONS
    """
}

process AGGREGATE_METRICS {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path metrics_files

    output:
    path "metrics.json", emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    aggregate_metrics.py \\
        ${metrics_files} \\
        --output metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        aggregate_metrics: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_JSON > metrics.json
    {
        "per_tool_metrics": [],
        "summary": {
            "weighted_f1": 0.0,
            "species_precision": 0.0,
            "species_recall": 0.0,
            "species_f1": 0.0,
            "objective_rank": "species",
            "aggregation": "best_of_n"
        }
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        aggregate_metrics: 1.0.0
    END_VERSIONS
    """
}

workflow BENCHMARKING {
    take:
    taxpasta_profiles   // channel: [ val(meta), path(taxpasta_tsv) ] - standardised profiles from TAXPASTA
    truth_profile       // path: gold standard profile in CAMI Bioboxes format
    taxdump_dir         // path: NCBI taxonomy dump directory
    tools_used          // val: comma-separated list of profiling tools used

    main:
    ch_versions = Channel.empty()

    // Step 1: Convert TAXPASTA output to CAMI Bioboxes profiling format
    CONVERT_TO_BIOBOXES(
        taxpasta_profiles,
        taxdump_dir,
    )
    ch_versions = ch_versions.mix(CONVERT_TO_BIOBOXES.out.versions.first())

    // Step 2: Evaluate each profile against gold standard
    EVALUATE_PROFILE(
        CONVERT_TO_BIOBOXES.out.profile,
        truth_profile,
        tools_used,
    )
    ch_versions = ch_versions.mix(EVALUATE_PROFILE.out.versions.first())

    // Step 3: Aggregate per-tool metrics into single summary for Stimulus
    ch_all_metrics = EVALUATE_PROFILE.out.metrics
        .map { _meta, metrics_file -> metrics_file }
        .collect()

    AGGREGATE_METRICS(ch_all_metrics)
    ch_versions = ch_versions.mix(AGGREGATE_METRICS.out.versions)

    emit:
    per_tool_metrics = EVALUATE_PROFILE.out.metrics      // channel: [ val(meta), path(metrics.json) ]
    metrics          = AGGREGATE_METRICS.out.metrics      // path: aggregated metrics.json
    versions         = ch_versions                        // channel: [ path(versions.yml) ]
}
