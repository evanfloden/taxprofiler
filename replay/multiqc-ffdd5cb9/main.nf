nextflow.enable.dsl = 2

include { MULTIQC } from '../../modules/nf-core/multiqc/main'

workflow {
    ch_multiqc_files = Channel
        .fromPath("${params.data}/*", checkIfExists: true)
        .filter { file -> !file.name.startsWith('.') }
        .collect()

    MULTIQC(
        ch_multiqc_files,
        file("${projectDir}/replay/multiqc-ffdd5cb9/multiqc_config.yml", checkIfExists: true),
        [],
        [],
        [],
        [],
    )
}
