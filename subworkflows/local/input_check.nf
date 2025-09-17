//
// Process pre-validated samplesheet data and get read channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // channel: pre-validated samplesheet data from nf-schema
    seq_center  // string: sequencing center for read group

    main:
    ch_versions = Channel.empty()

    samplesheet
        .map { create_fastq_channel(it, seq_center) }
        .set { reads }

    emit:
    reads                 // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(data_list, String seq_center) {
    // nf-schema returns: [meta, fastq_1, fastq_2, replicate, antibody, control, control_replicate]
    // Extract the parts
    def meta_from_schema = data_list[0]
    def fastq_1 = data_list[1]
    def fastq_2 = data_list[2]
    def replicate = data_list[3] ?: 1
    def antibody = data_list[4] ?: ''
    def control = data_list[5] ?: ''
    def control_replicate = data_list[6] ?: 1

    def meta = [:]
    meta.id         = meta_from_schema.id
    meta.single_end = (fastq_2 == null || fastq_2.toString().trim() == '' || fastq_2 instanceof List && fastq_2.isEmpty())
    meta.antibody   = antibody
    meta.control    = control
    meta.replicate  = replicate
    meta.control_replicate = control_replicate

    def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = read_group

    // File paths are already validated by nf-schema
    def fastq_meta = []
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(fastq_1) ] ]
    } else {
        fastq_meta = [ meta, [ file(fastq_1), file(fastq_2) ] ]
    }
    return fastq_meta
}