#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Thread-safe list to store completion messages
def completion_log = [].asSynchronized()

// Define parameters with default values
params.pod5_dir = false
params.ref_genome = false
params.ref_transcriptome = false

params.outdir = './results'
params.help = false
params.model = 'sup'   // Default to 'sup' (super accuracy) model
params.use_gpu = 'cuda:all'     // GPU specification: 'cuda:all' for all GPUs, or specify GPU IDs like 'cuda:0', 'cuda:0,1', etc.
params.cpus = 4          // Default number of CPUs to use

// Modification Calling Parameters
params.mods = false //"m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG"  // Specify modifications to call, e.g., 'm6A 5mC pseU'.

// Modkit Pileup Parameters
params.modkit_filter_threshold = 0.8
params.modkit_mod_threshold    = 0.98

// Modkit mod name normalization and threshold flag builder
def MODKIT_MOD_MAP = [
    'm6a':   'a',
    '6ma':   'a',
    'am':    '69426',
    '2omea': '69426',
    'inosine': '17596',
    'ino':     '17596',
    'pseu':    '17802',
    'pseudouridine': '17802',
    'um':    '19227',
    '2omeu': '19227',
    'gm':    '19229',
    '2omeg': '19229',
    'cm':    '19228',
    '2omec': '19228',
    'm5c':   'm',
    '5mc':   'm'
]

def normalizeModName(String mod) {
    if (!mod) {
        return null
    }
    return mod.toLowerCase()
              .replaceAll(/\s+/, '')
              .replaceAll(/[^a-z0-9]/, '')
}

def buildModThresholdFlags(def mods, Map modMap, def threshold) {
    if (!mods || mods == "false") {
        return ""
    }
    def tokens = mods.toString()
                     .split(',') as List
    tokens = tokens.collectMany { (it.split('_') as List) }
                     .collect { normalizeModName(it) }
                     .findAll { it }

    def flags = tokens.collect { key ->
        def code = modMap[key]
        code ? "--mod-threshold ${code}:${threshold}" : null
    }.findAll { it }

    return flags.unique().join(' ')
}

def buildModifiedBasesFlag(def mods) {
    if (!mods || mods == "false") {
        return ""
    }
    def tokens = mods.toString()
                     .split(',') as List
    tokens = tokens.collectMany { (it.split('_') as List) }
                     .collect { it?.trim() }
                     .findAll { it }
    return tokens ? "--modified-bases ${tokens.join(' ')}" : ""
}

// Help message
if (params.help) {
    log.info """
    Nanopore RNA Modification Analysis Pipeline
    =========================================
    This pipeline performs basecalling with integrated modification calling, alignment,
    and modification pileup analysis using Dorado and Modkit.
    
    Usage:
    nextflow run main.nf --pod5_dir <pod5_directory> --ref_genome <reference.fasta> [options]
    
    Required Arguments:
      --pod5_dir      Path to the directory containing POD5 files (use quotes and wildcards for multiple samples).
      --ref_genome    Path to the reference genome in FASTA format.
      --ref_transcriptome    Path to the reference transcriptome in FASTA format.
    
    Optional Arguments:
      --outdir              Base path for output directories. Results will be organized by file type (base, featureC, salmon, etc.). Default: ./results
      --model               Basecalling model accuracy. Default: 'sup'
      --use_gpu             GPU specification for basecalling. Default: 'cuda:all'
                            - 'cuda:all': Use all available GPUs
                            - 'cuda:0': Use GPU 0 only
                            - 'cuda:0,1': Use GPUs 0 and 1
                            - Specify any valid GPU ID(s) separated by commas
      --cpus                Number of CPUs to use. Default: 4
      
    Modification Calling Options:
      --mods                Space-separated list of modifications to call (e.g., 'm6A').
                            Set to false to disable modification calling. Default: false
      --modkit_threshold    Probability threshold for Modkit pileup. Default: 0.98
      --motif_str           Motif and offset for Modkit analysis (e.g., 'DRACH 2'). Default: false (all-context)

      --help                Display this help message.
    """
    exit 0
}

// Check for required parameters
if (!params.pod5_dir) {
    exit 1, "Input POD5 directory not specified! Please use --pod5_dir"
}
if (!params.ref_genome) {
    exit 1, "Reference genome not specified! Please use --ref_genome"
}
if (!params.ref_transcriptome) {
    exit 1, "Reference transcriptome not specified! Please use --ref_transcriptome"}



// ---
// Processes
// ---


// Basecalling with Dorado, including modification calling if specified
process DORADO_BASECALL {
    publishDir "${params.outdir}/base", mode: 'copy'
    
    input:
    tuple val(meta), path(pod5_dir)
    
    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    
    script:
    def gpu_flag = "--device ${params.use_gpu}"
    def mod_flag = (params.mods && params.mods != "false") ? "--modified-bases ${params.mods.replaceAll(',', ' ')}" : ''

    """
    dorado basecaller \\
        ${gpu_flag} \\
        "${params.model}" \\
        "${pod5_dir}" \\
        ${mod_flag} -n 10000 > "${meta.id}.bam"
    """
}

// Genome alignment with Minimap2 conserving ML and MM tags for modifications
process MINIMAP2_ALIGN_GENOME {
    publishDir "${params.outdir}/base", mode: 'copy'
    
    input:
    tuple val(meta), path(bam_file)
    path ref_genome
    
    output:
    tuple val(meta), path("${meta.id}.aligned.sorted.bam"), path("${meta.id}.aligned.sorted.bam.bai"), emit: bam_bai
    
    script:
    // Use 'splice-ont' for genome alignment of RNA reads
    """
    # Extract @RG lines from original Dorado BAM
    samtools view -H "${bam_file}" | grep "^@RG" > original_rg.txt

    # This preserves modification tags (mm, ml) essential for modification analysis
    samtools fastq -T "*" "${bam_file}" | minimap2 -y --MD -ax splice-ont -t ${task.cpus} "${ref_genome}" - | \\
    samtools view -bh -F 260 | \\
    samtools sort -@ ${task.cpus} -o temp.bam

    # Merge headers: aligned BAM header + original @RG lines
    samtools view -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam

    # Apply new header
    samtools reheader new_header.sam temp.bam > "${meta.id}.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    
    # Index the BAM file
    samtools index "${meta.id}.aligned.sorted.bam"
    """
}

// Transcriptome alignment with Minimap2 conserving ML and MM tags for modifications
process MINIMAP2_ALIGN_TRANSCRIPTOME {
    publishDir "${params.outdir}/base", mode: 'copy'
    
    input:
    tuple val(meta), path(bam_file)
    path ref_transcriptome
    
    output:
    tuple val(meta), path("${meta.id}.transcriptome.aligned.sorted.bam"), path("${meta.id}.transcriptome.aligned.sorted.bam.bai"), emit: bam_bai
    
    script:
    // Use 'map-ont' preset for transcriptome alignment of RNA reads
    """
    # Extract @RG lines from original Dorado BAM
    samtools view -H "${bam_file}" | grep "^@RG" > original_rg.txt

    samtools fastq -T "*" "${bam_file}" | minimap2 -y --MD -ax map-ont -t ${task.cpus} "${ref_transcriptome}" - | \\
    samtools view -bh -F 260 | \\
    samtools sort -@ ${task.cpus} -o temp.bam

    # Merge headers: aligned BAM header + original @RG lines
    samtools view -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam

    # Apply new header
    samtools reheader new_header.sam temp.bam > "${meta.id}.transcriptome.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    
    # Index the BAM file
    samtools index "${meta.id}.transcriptome.aligned.sorted.bam"
    """
}

//NanoComp process for QC check with NanoComp Stats
process NANOCOMP {
    publishDir "${params.outdir}/nanocomp", mode: 'copy'

    
    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bai)
    tuple val(meta2), path(genome_bam), path(genome_bai)

    output:
    tuple val(meta), path("nanocomp_report_${meta.id}"), emit: txt

    script:
    """
    # Run NanoComp comparing transcriptome and genome alignments for the same sample
    NanoComp \
        --bam ${transcriptome_bam} ${genome_bam} \\
        --plot false \\
        --outdir "nanocomp_report_${meta.id}"
    """
}

// Modkit pileup for genome-aligned BAM
process MODKIT_PILEUP_GENOME {
    publishDir "${params.outdir}/modkitGenome", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path ref_genome
    
    output:
    tuple val(meta), path("${meta.id}.genome.bed"), emit: bed
    tuple val(meta), path("${meta.id}.genome.log"), emit: log
    
    script:
    // Parse modifications called during basecalling and generate mod-threshold flags
    def mod_threshold_flags = buildModThresholdFlags(params.mods, MODKIT_MOD_MAP, params.modkit_mod_threshold)
    def modified_bases_flag = buildModifiedBasesFlag(params.mods)
    
    """
    modkit pileup \\
        ${bam} \\
        ${meta.id}.genome.bed \\
        --ref ${ref_genome} \\
        --threads ${task.cpus} \\
        --filter-threshold ${params.modkit_filter_threshold} \\
        ${modified_bases_flag} \
        ${mod_threshold_flags} \\
        --log-filepath "${meta.id}.genome.log" \\
        --bedrmod
    """
}

// Modkit pileup for transcriptome-aligned BAM
process MODKIT_PILEUP_TRANSCRIPTOME {
    publishDir "${params.outdir}/modkitTranscriptome", mode: 'copy', pattern: "*.{bed,log}"
    
    input:
    tuple val(meta), path(bam), path(bai)
    path ref_transcriptome
    
    output:
    tuple val(meta), path("${meta.id}.transcriptome.bed"), emit: bed
    tuple val(meta), path("${meta.id}.transcriptome.log"), emit: log
    
    script:
    // Parse modifications called during basecalling and generate mod-threshold flags
    def mod_threshold_flags = buildModThresholdFlags(params.mods, MODKIT_MOD_MAP, params.modkit_mod_threshold)
    def modified_bases_flag = buildModifiedBasesFlag(params.mods)
    
    """
    modkit pileup \\
        ${bam} \\
        ${meta.id}.transcriptome.bed \\
        --ref ${ref_transcriptome} \\
        --threads ${task.cpus} \\
        --filter-threshold ${params.modkit_filter_threshold} \\
        ${modified_bases_flag} \
        ${mod_threshold_flags} \\
        --log-filepath "${meta.id}.transcriptome.log" \\
        --bedrmod
    """
}



// ---
// Workflow
// ---

workflow {
    // Create input channels with metadata
    // Fix the input channel to properly handle directories of POD5 files
    ch_pod5_input = Channel.fromPath(params.pod5_dir, type: 'dir')
                        .map { dir -> tuple([id: dir.getBaseName()], dir) }
    
    ch_ref_genome = file(params.ref_genome)
    ch_ref_transcriptome = file(params.ref_transcriptome)

    // Print debug information
    ch_pod5_input.view { meta, dir -> "Processing POD5 directory: ${dir} with ID: ${meta.id}" }

    // 1. Basecalling with optional modification calling
    DORADO_BASECALL(ch_pod5_input)
    DORADO_BASECALL.out.bam.subscribe { meta, file -> 
        def msg = "Completed Dorado Basecalling for: ${meta.id}\n  -> Output: ${params.outdir}/base/${file.name}"
        log.info msg
    }

    // 2. Genome alignment (always required for featureCounts)
    MINIMAP2_ALIGN_GENOME(DORADO_BASECALL.out.bam, ch_ref_genome)
    MINIMAP2_ALIGN_GENOME.out.bam_bai.subscribe { meta, bam, bai -> 
        def msg = "Completed Genome Alignment for: ${meta.id}\n  -> Output: ${params.outdir}/base/${bam.name}"
        log.info msg
    }

    // 3. Transcriptome alignment
    MINIMAP2_ALIGN_TRANSCRIPTOME(DORADO_BASECALL.out.bam, ch_ref_transcriptome)
    MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam_bai.subscribe { meta, bam, bai -> 
        def msg = "Completed Transcriptome Alignment for: ${meta.id}\n  -> Output: ${params.outdir}/base/${bam.name}"
        log.info msg
    }


    // 7. Modification pileup (only if mods were called)

    MODKIT_PILEUP_GENOME(MINIMAP2_ALIGN_GENOME.out.bam_bai, ch_ref_genome)
    // Run transcriptome acceleration before regular modkit pileup
    MODKIT_PILEUP_TRANSCRIPTOME(MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam_bai, ch_ref_transcriptome)

    
    // Log completion for modification steps
    MODKIT_PILEUP_GENOME.out.bed.subscribe { meta, file -> 
        def msg = "Completed Modkit Genome Pileup for: ${meta.id}\n  -> Output: ${params.outdir}/modkitGenome/${file.name}"
        log.info msg
        completion_log.add(msg)
    }
    MODKIT_PILEUP_TRANSCRIPTOME.out.bed.subscribe { meta, file -> 
        def msg = "Completed Modkit Transcriptome Pileup for: ${meta.id}\n  -> Output: ${params.outdir}/modkitTranscriptome/${file.name}"
        log.info msg
        completion_log.add(msg)
    }
}

workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """

    log.info "\nGenerated Output Files:"
    log.info "======================="
    completion_log.each { msg -> log.info msg }
}