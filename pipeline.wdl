version 1.0

workflow nanopore_modcall {
    input {
        Array[File] pod5_files
        File ref_genome
        File ref_transcriptome
        String use_gpu = "cuda:all"
        Int cpus = 64
    }

    scatter (pod5_file in pod5_files) {
        String sample_id = basename(pod5_file, ".pod5")

        call DoradoBasecall {
            input:
                pod5_file = pod5_file,
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
    }

    output {
        Array[File] basecall_bams = DoradoBasecall.bam
        Array[File] genome_bams = MinimapGenome.aligned_bam
        Array[File] genome_bais = MinimapGenome.bai
        Array[File] transcriptome_bams = MinimapTranscriptome.aligned_bam
        Array[File] transcriptome_bais = MinimapTranscriptome.bai
        Array[File] nanocomp_reports = NanoCompQC.report_tar
        Array[File] genome_beds = ModkitPileupGenome.bed
        Array[File] genome_logs = ModkitPileupGenome.log
        Array[File] transcriptome_beds = ModkitPileupTranscriptome.bed
        Array[File] transcriptome_logs = ModkitPileupTranscriptome.log
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
        --device ~{use_gpu} \
        "sup,inosine_m6A_2OmeA,m5C_2OmeC,pseU_2OmeU,2OmeG" \
        "~{pod5_file}" \
        -n 10000 > "~{sample_id}.bam"
    >>>

    output {
        File bam = "~{sample_id}.bam"
    }

    runtime {
        cpu: cpus
        # Use digest form with @ to avoid manifest lookup failures
        docker: "ontresearch/dorado@sha256:6156fdbb48ff13fbb141b4f1fc6e6300f05221ecfdd114fc8323a0e38296c1dd"
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

    samtools view -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq -T "*" "~{bam}" | minimap2 -y --MD -ax splice-ont -t ~{cpus} "~{ref_genome}" - | \
    samtools view -bh -F 260 | \
    samtools sort -@ ~{cpus} -o temp.bam

    samtools view -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "~{sample_id}.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index "~{sample_id}.aligned.sorted.bam"
    >>>

    output {
        File aligned_bam = "~{sample_id}.aligned.sorted.bam"
        File bai = "~{sample_id}.aligned.sorted.bam.bai"
    }

    runtime {
        cpu: cpus
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

    samtools view -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq -T "*" "~{bam}" | minimap2 -y --MD -ax map-ont -t ~{cpus} "~{ref_transcriptome}" - | \
    samtools view -bh -F 260 | \
    samtools sort -@ ~{cpus} -o temp.bam

    samtools view -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "~{sample_id}.transcriptome.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index "~{sample_id}.transcriptome.aligned.sorted.bam"
    >>>

    output {
        File aligned_bam = "~{sample_id}.transcriptome.aligned.sorted.bam"
        File bai = "~{sample_id}.transcriptome.aligned.sorted.bam.bai"
    }

    runtime {
        cpu: cpus
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
             --outdir "nanocomp_report_~{sample_id}" \
             --logfile "nanocomp_report_~{sample_id}/nanocomp.log" \
             --report "nanocomp_report_~{sample_id}/nanocomp.html"

    tar -czf "nanocomp_report_~{sample_id}.tar.gz" "nanocomp_report_~{sample_id}"
    >>>

    output {
        File report_tar = "nanocomp_report_~{sample_id}.tar.gz"
    }

    runtime {
        cpu: cpus
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
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
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
        --bedrmod
    >>>

    output {
        File bed = "~{sample_id}.transcriptome.bed"
        File log = "~{sample_id}.transcriptome.log"
    }

    runtime {
        cpu: cpus
        docker: "ontresearch/modkit:sha489d708a48c66368e5d1e118538e5dca68203a64"
    }
}
