/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BENCHMARKING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Compares taxonomic profiles against a gold standard using OPAL.
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

process BENCHMARK_OPAL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::cami-opal=1.0.14"
    container "community.wave.seqera.io/library/cami-opal:23aa9c620f30c4b3"

    input:
    tuple val(meta), path(query_profile)
    path truth_profile

    output:
    tuple val(meta), path("opal_results/"), emit: results
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p opal_results

    opal.py \\
        -g ${truth_profile} \\
        ${query_profile} \\
        -l "${meta.id}" \\
        -o opal_results/ \\
    || true

    # If OPAL failed, create a minimal results.tsv with zero metrics
    if [ ! -f opal_results/results.tsv ]; then
        mkdir -p opal_results
        printf "tool\\trank\\tmetric\\tvalue\\n" > opal_results/results.tsv
        printf "${meta.id}\\tspecies\\tpurity (precision)\\t0.0\\n" >> opal_results/results.tsv
        printf "${meta.id}\\tspecies\\tcompleteness (recall)\\t0.0\\n" >> opal_results/results.tsv
        printf "${meta.id}\\tspecies\\tF1 score\\t0.0\\n" >> opal_results/results.tsv
        printf "${meta.id}\\tspecies\\tL1 norm error\\t2.0\\n" >> opal_results/results.tsv
        printf "${meta.id}\\tspecies\\tWeighted UniFrac error\\t16.0\\n" >> opal_results/results.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opal: \$(opal.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo '1.0.14')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p opal_results
    printf "tool\\trank\\tmetric\\tvalue\\n" > opal_results/results.tsv
    printf "${meta.id}\\tspecies\\tpurity (precision)\\t0.85\\n" >> opal_results/results.tsv
    printf "${meta.id}\\tspecies\\tcompleteness (recall)\\t0.80\\n" >> opal_results/results.tsv
    printf "${meta.id}\\tspecies\\tF1 score\\t0.824\\n" >> opal_results/results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opal: 1.0.14
    END_VERSIONS
    """
}

process EXTRACT_METRICS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(opal_results)
    val tools_used

    output:
    tuple val(meta), path("metrics.json"), emit: metrics
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_metrics.py \\
        --opal-results ${opal_results} \\
        --sample ${meta.id} \\
        --tools "${tools_used}" \\
        --output metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        extract_metrics: 1.0.0
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
        extract_metrics: 1.0.0
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

    // Step 2: Run OPAL evaluation against gold standard
    BENCHMARK_OPAL(
        CONVERT_TO_BIOBOXES.out.profile,
        truth_profile,
    )
    ch_versions = ch_versions.mix(BENCHMARK_OPAL.out.versions.first())

    // Step 3: Extract metrics from OPAL output
    EXTRACT_METRICS(
        BENCHMARK_OPAL.out.results,
        tools_used,
    )
    ch_versions = ch_versions.mix(EXTRACT_METRICS.out.versions.first())

    emit:
    metrics  = EXTRACT_METRICS.out.metrics   // channel: [ val(meta), path(metrics.json) ]
    versions = ch_versions                    // channel: [ path(versions.yml) ]
}
