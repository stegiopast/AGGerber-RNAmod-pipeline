version 1.0

workflow ont_mRNA_pilot {
    input {
        Array[File] pod5_files
        String sample_id
        File ref_genome
        File ref_transcriptome
        String use_gpu = "cuda:all"
        Int cpus = 64
    }

    if (length(pod5_files) > 1) {
        call Pod5Merge {
            input:
                pod5_files = pod5_files,
                sample_id = sample_id,
                cpus = cpus
        }
    }

    File pod5_input = select_first([Pod5Merge.merged_pod5, pod5_files[0]])

    call DoradoBasecall {
        input:
            pod5_file = pod5_input,
            sample_id = sample_id,
            use_gpu = use_gpu,
            cpus = cpus
    }

    call MinimapGenome {
        input:
            bam = DoradoBasecall.bam,
            sample_id = sample_id,
            ref_genome = ref_genome,
            cpus = cpus
    }

    call MinimapTranscriptome {
        input:
            bam = DoradoBasecall.bam,
            sample_id = sample_id,
            ref_transcriptome = ref_transcriptome,
            cpus = cpus
    }

    call NanoCompQC {
        input:
            genome_bam = MinimapGenome.aligned_bam,
            transcriptome_bam = MinimapTranscriptome.aligned_bam,
            sample_id = sample_id,
            cpus = cpus
    }

    call ModkitPileupGenome {
        input:
            bam = MinimapGenome.aligned_bam,
            bai = MinimapGenome.bai,
            sample_id = sample_id,
            ref_genome = ref_genome,
            cpus = cpus
    }

    call ModkitPileupTranscriptome {
        input:
            bam = MinimapTranscriptome.aligned_bam,
            bai = MinimapTranscriptome.bai,
            sample_id = sample_id,
            ref_transcriptome = ref_transcriptome,
            cpus = cpus
    }

    output {
        File basecall_bam = DoradoBasecall.bam
        File genome_bam = MinimapGenome.aligned_bam
        File genome_bai = MinimapGenome.bai
        File transcriptome_bam = MinimapTranscriptome.aligned_bam
        File transcriptome_bai = MinimapTranscriptome.bai
        File nanocomp_report = NanoCompQC.report_tar
        File genome_bed = ModkitPileupGenome.bed
        File genome_log = ModkitPileupGenome.log
        File transcriptome_bed = ModkitPileupTranscriptome.bed
        File transcriptome_log = ModkitPileupTranscriptome.log
    }
}

task Pod5Merge {
    input {
        Array[File] pod5_files
        String sample_id
        Int cpus
    }

    command <<<!
    set -euo pipefail

    pod5 merge -t ~{cpus} -D ~{sep=' ' pod5_files} -o "~{sample_id}.merged.pod5"
    >>>

    output {
        File merged_pod5 = "~{sample_id}.merged.pod5"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds) Maximum Allocation time = (~48h = 172800)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        # Use digest form with @ to avoid manifest lookup failures
        docker: "ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
    }
}

task DoradoBasecall {
    input {
        File pod5_file
        String sample_id
        String use_gpu
        Int cpus
    }

    command <<<
    set -euo pipefail

    dorado basecaller \
        --device "~{use_gpu}" \
        "sup,inosine_m6A_2OmeA,m5C_2OmeC,pseU_2OmeU,2OmeG" \
        --estimate-poly-a \
        --emit-moves \
        "~{pod5_file}" > "~{sample_id}.bam"
    >>>

    output {
        # When --emit-moves is used, Dorado embeds the move table in the BAM.
        File bam = "~{sample_id}.bam"
    }

    runtime {
        cpu: cpus
        gpu: true
        gpuCount: 1
        memory: "512GB"
        maxRunTime: 172800 #48 hours (48 * 3600 seconds)
        runtime_minutes: 2820 #47 hours (47 * 60 minutes)
        # Use digest form with @ to avoid manifest lookup failures
        docker: "ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
    }
}

task MinimapGenome {
    input {
        File bam
        String sample_id
        File ref_genome
        Int cpus
    }

    command <<<
    set -euo pipefail

    samtools view --threads ~{cpus} -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq --threads ~{cpus} -T "*" "~{bam}" | minimap2 -y --MD -ax splice -uf -k14 -t ~{cpus} "~{ref_genome}" - | \
    samtools sort --threads ~{cpus} | \
    samtools view --threads ~{cpus} -bh -F 260 -o temp.bam

    samtools view --threads ~{cpus} -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "~{sample_id}.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index -@ ~{cpus} "~{sample_id}.aligned.sorted.bam"
    >>>

    output {
        File aligned_bam = "~{sample_id}.aligned.sorted.bam"
        File bai = "~{sample_id}.aligned.sorted.bam.bai"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        # Match the working Nextflow pipeline tag
        docker: "nanozoo/minimap2:2.28--9e3bd01"
    }
}

task MinimapTranscriptome {
    input {
        File bam
        String sample_id
        File ref_transcriptome
        Int cpus
    }

    command <<<
    set -euo pipefail

    samtools view --threads ~{cpus} -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq --threads ~{cpus} -T "*" "~{bam}" | minimap2 -y --MD -ax map-ont -t ~{cpus} "~{ref_transcriptome}" - | \
    samtools sort --threads ~{cpus} | \
    samtools view --threads ~{cpus} -bh -F 260 -o temp.bam

    samtools view --threads ~{cpus} -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "~{sample_id}.transcriptome.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index -@ ~{cpus} "~{sample_id}.transcriptome.aligned.sorted.bam"
    >>>

    output {
        File aligned_bam = "~{sample_id}.transcriptome.aligned.sorted.bam"
        File bai = "~{sample_id}.transcriptome.aligned.sorted.bam.bai"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        # Match the working Nextflow pipeline tag
        docker: "nanozoo/minimap2:2.28--9e3bd01"
    }
}

task NanoCompQC {
    input {
        File genome_bam
        File transcriptome_bam
        String sample_id
        Int cpus
    }

    command <<<
    set -euo pipefail

    mkdir -p "nanocomp_report_~{sample_id}"
    NanoComp --bam "~{genome_bam}" "~{transcriptome_bam}" \
    --names "genome" "transcriptome" \
    --threads ~{cpus} \
    --outdir "nanocomp_report_~{sample_id}"

    tar -czf "nanocomp_report_~{sample_id}.tar.gz" "nanocomp_report_~{sample_id}"
    >>>

    output {
        File report_tar = "nanocomp_report_~{sample_id}.tar.gz"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 43200 #12 hours (12 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        docker: "luxendr13/nanocomp:0.6.0"
    }
}

task ModkitPileupGenome {
    input {
        File bam
        File bai
        String sample_id
        File ref_genome
        Int cpus
    }

    command <<<
    set -euo pipefail

    modkit pileup \
        ~{bam} \
        ~{sample_id}.genome.bed \
        --ref ~{ref_genome} \
        --threads ~{cpus} \
        --filter-threshold 0.8 \
        --modified-bases m5C 2OmeC inosine m6A 2OmeA pseU 2OmeU 2OmeG \
        --mod-threshold m:0.98 \
        --mod-threshold 19228:0.98 \
        --mod-threshold 17596:0.98 \
        --mod-threshold a:0.98 \
        --mod-threshold 69426:0.98 \
        --mod-threshold 17802:0.98 \
        --mod-threshold 19227:0.98 \
        --mod-threshold 19229:0.98 \
        --log-filepath "~{sample_id}.genome.log" \
        --bedrmod
    >>>

    output {
        File bed = "~{sample_id}.genome.bed"
        File log = "~{sample_id}.genome.log"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
        failOnStderr: false
    }
}

task ModkitPileupTranscriptome {
    input {
        File bam
        File bai
        String sample_id
        File ref_transcriptome
        Int cpus
    }

    command <<<
    set -euo pipefail

    modkit pileup \
        ~{bam} \
        ~{sample_id}.transcriptome.bed \
        --ref ~{ref_transcriptome} \
        --threads ~{cpus} \
        --filter-threshold 0.8 \
        --modified-bases m5C 2OmeC inosine m6A 2OmeA pseU 2OmeU 2OmeG \
        --mod-threshold m:0.98 \
        --mod-threshold 19228:0.98 \
        --mod-threshold 17596:0.98 \
        --mod-threshold a:0.98 \
        --mod-threshold 69426:0.98 \
        --mod-threshold 17802:0.98 \
        --mod-threshold 19227:0.98 \
        --mod-threshold 19229:0.98 \
        --log-filepath "~{sample_id}.transcriptome.log" \
        --preload-references \
        --bedrmod
    >>>

    output {
        File bed = "~{sample_id}.transcriptome.bed"
        File log = "~{sample_id}.transcriptome.log"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
        failOnStderr: false
    }
}
