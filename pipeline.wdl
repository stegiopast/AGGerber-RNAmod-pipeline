version 1.0

workflow ont_mRNA_pilot {
    input {
        Array[File] pod5_files
        String sample_id
        File ref_genome
        File ref_transcriptome
        String use_gpu
        Int cpus
        Int gpus
        Boolean barcoded = false
        String barcoding_model = "b96_RNA004"
        String basecalling_type = "sup"
    }

    # if (length(pod5_files) > 1) {
    #     call Pod5Merge {
    #         input:
    #             pod5_files = pod5_files,
    #             sample_id = sample_id,
    #             cpus = cpus
    #     }
    # }

    # File pod5_input = select_first([Pod5Merge.merged_pod5, pod5_files[0]])

    # call DoradoBasecall {
    #     input:
    #         pod5_file = pod5_input,
    #         sample_id = sample_id,
    #         use_gpu = use_gpu,
    #         cpus = cpus
    # }

    call Seqtagger {
        input:
            pod5_files = pod5_files,
            sample_id = sample_id,
            use_gpu = use_gpu,
            cpus = cpus,
            gpus = gpus,
            barcoded = barcoded,
            barcoding_model = barcoding_model
    }
    
    call DoradoBasecall {
        input:
            barcode_table = Seqtagger.barcode_table,
            pod5_files = pod5_files,
            sample_id = sample_id,
            use_gpu = use_gpu,
            cpus = cpus,
            gpus = gpus,
            basecalling_type = basecalling_type
        }

    scatter (chunk in DoradoBasecall.bam) {

        call MinimapGenome {
            input:
                bam = chunk,
                sample_id = sample_id,
                ref_genome = ref_genome,
                cpus = cpus
        }

        call MinimapTranscriptome {
            input:
                bam = chunk,
                sample_id = sample_id,
                ref_transcriptome = ref_transcriptome,
                cpus = cpus
        }

        if (MinimapGenome.aligned_count > 0) {
        call NanoCompQC {
            input:
                genome_bam = MinimapGenome.aligned_bam,
                transcriptome_bam = MinimapTranscriptome.aligned_bam,
                sample_id = sample_id,
                cpus = cpus
            }
        }

        if (basecalling_type == "sup" && MinimapGenome.aligned_count > 0) { 
            call ModkitPileupGenome {
            input:
                bam = MinimapGenome.aligned_bam,
                bai = MinimapGenome.aligned_bai,
                sample_id = sample_id,
                ref_genome = ref_genome,
                cpus = cpus
            }
        }    

        if (basecalling_type == "sup" && MinimapTranscriptome.aligned_count > 0) {
            call ModkitPileupTranscriptome {
            input:
                bam = MinimapTranscriptome.aligned_bam,
                bai = MinimapTranscriptome.aligned_bai,
                sample_id = sample_id,
                ref_transcriptome = ref_transcriptome,
                cpus = cpus
            }
        } 

    }
    output {
        File barcode_table = Seqtagger.barcode_table
        Array[File] barcode_pdfs = Seqtagger.pdfs
        Array[File] basecall_bam = DoradoBasecall.bam
        Array[File] genome_bam = MinimapGenome.aligned_bam
        Array[File] genome_bai = MinimapGenome.aligned_bai
        Array[File] transcriptome_bam = MinimapTranscriptome.aligned_bam
        Array[File] transcriptome_bai = MinimapTranscriptome.aligned_bai
        Array[File?] nanocomp_report = NanoCompQC.report_tar
        Array[File?] genome_bed = ModkitPileupGenome.bed
        Array[File?] genome_log = ModkitPileupGenome.log
        Array[File?] transcriptome_bed = ModkitPileupTranscriptome.bed
        Array[File?] transcriptome_log = ModkitPileupTranscriptome.log
    }
    }


# task Pod5Merge {
#     input {
#         Array[File] pod5_files
#         String sample_id
#         Int cpus
#     }

#     command <<<!
#     set -euo pipefail

#     pod5 merge -t ~{cpus} -D ~{sep=' ' pod5_files} -o "~{sample_id}.merged.pod5"
#     >>>

#     output {
#         File merged_pod5 = "~{sample_id}.merged.pod5"
#     }

#     runtime {
#         cpu: cpus
#         memory: "128GB"
#         maxRunTime: 86400 #24 hours (24 * 3600 seconds) Maximum Allocation time = (~48h = 172800)
#         runtime_minutes: 1440 #24 hours (24 * 60 minutes)
#         # Use digest form with @ to avoid manifest lookup failures
#         docker: "ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
#     }
# }

task Seqtagger {
    input {
        Array[File] pod5_files
        String sample_id
        String use_gpu
        Int cpus
        Int gpus
        Boolean barcoded
        String barcoding_model
    }

    command <<<
    set -euo pipefail

    if [ ~{barcoded} = true ]; then
    for pod5 in ~{sep=' ' pod5_files}; do
        filename=$(basename "$pod5" .pod5)
        mkdir -p ./$filename
        mRNA -k /opt/app/models/~{barcoding_model} -r -i $pod5 -o ./$filename/
        zcat ./$filename/*.demux.tsv.gz | awk 'BEGIN { FS = OFS = "\t" } $5 >= 50 { print }' > ./$filename/$filename.tsv
        mv ./$filename/$filename.tsv ./
        mv ./$filename/*.pdf ./
        rm -rf ./$filename
        head -n 1 ./$filename.tsv > header.tsv
    done
    tail -n +2 -q ./*.tsv | sort -k1,1 -k2,2n | cat header.tsv - > "demux_output.tsv"
    else
        echo "Skip" > "demux_output.tsv"
        echo "Skip" > "no_barcoding.pdf"
    fi
    >>>

    output {
        File barcode_table = "demux_output.tsv"
        Array[File] pdfs = glob("*.pdf")
    }

    runtime {
        cpu: cpus
        gpu: true
        gpuCount: gpus
        memory: "128GB"
        maxRunTime: 86400 #24 hours (24 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        docker: "lpryszcz/seqtagger:latest"
    }
}

task DoradoBasecall {
    input {
        File barcode_table
        Array[File] pod5_files
        String sample_id
        String use_gpu
        Int cpus
        Int gpus
        String basecalling_type
    }

    command <<<
    set -euo pipefail

    mkdir -p bam_temp
    mkdir -p final_bams
    for pod5 in ~{sep=' ' pod5_files}; do
        filename=$(basename "$pod5" .pod5)
        if [ ~{basecalling_type} == "sup" ]; then
        dorado basecaller \
        --device ~{use_gpu} \
        "sup,inosine_m6A_2OmeA,m5C_2OmeC,pseU_2OmeU,2OmeG" \
        --estimate-poly-a \
        --emit-moves $pod5 > "bam_temp/${filename}.bam"
        else
        dorado basecaller \
        --device ~{use_gpu} \
        "fast" \
        --estimate-poly-a \
        --emit-moves $pod5 > "bam_temp/${filename}.bam"
        fi
    done
    ls bam_temp/*.bam > bam_list.txt
    samtools merge -@ ~{cpus} -o "~{sample_id}.bam" -b bam_list.txt
    samtools index -b -@ ~{cpus} "~{sample_id}.bam"
    if grep -q "^Skip$" ~{barcode_table}; then
        echo "No barcoding performed, skipping barcode table cleanup."
        mv "~{sample_id}.bam" "final_bams/"
        mv "~{sample_id}.bam.bai" "final_bams/"
    else
        mkdir -p "barcode_tables"
        mkdir -p "final_bams"
        awk -F'\t' -v outdir="barcode_tables" 'NR > 1 {num=sprintf("%02d",$3); print > (outdir "/barcode" num ".tsv")}' ~{barcode_table}
        for tsv in barcode_tables/*.tsv; do
            filename=$(basename "$tsv" .tsv)
            cut -f 1 "$tsv" > "barcode_tables/${filename}.read_ids.txt"
            samtools view -@ ~{cpus} -bh -N "barcode_tables/${filename}.read_ids.txt" "~{sample_id}.bam" | samtools sort -@ ~{cpus} -o "${filename}.bam"
            samtools index -b -@ ~{cpus} "${filename}.bam"
            mv "${filename}.bam" "final_bams/"
            mv "${filename}.bam.bai" "final_bams/"
        done
        rm -rf "barcode_tables"
        rm "~{sample_id}.bam"
    fi
    rm bam_list.txt
    rm -rf "bam_temp"
    >>>

    output {
        # When --emit-moves is used, Dorado embeds the move table in the BAM.
        Array[File] bam = glob("final_bams/*.bam")
        Array[File] bai = glob("final_bams/*.bam.bai")
    }

    runtime {
        cpu: cpus
        gpu: true
        gpuCount: gpus
        memory: "512GB"
        maxRunTime: 172800 #48 hours (48 * 3600 seconds)
        runtime_minutes: 2820 #47 hours (47 * 60 minutes)
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

    filename=$(basename "~{bam}" .bam)
    samtools view --threads ~{cpus} -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq --threads ~{cpus} -T "*" "~{bam}" | minimap2 -y --MD -ax splice -uf -k14 -t ~{cpus} "~{ref_genome}" - | \
    samtools sort --threads ~{cpus} | \
    samtools view --threads ~{cpus} -bh -F 260 -o temp.bam

    samtools view --threads ~{cpus} -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "${filename}.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index -b -@ ~{cpus} "${filename}.aligned.sorted.bam"
    if [ -s "${filename}.aligned.sorted.bam" ]; then
        echo "Alignment successful for ${filename}.bam"
        samtools view --threads ~{cpus} -c "${filename}.aligned.sorted.bam" > aligned_count.txt
    else
        echo "Alignment failed for ${filename}.bam, creating empty BAM file."
        echo 0 > aligned_count.txt
        samtools view -hb -o "${filename}.aligned.sorted.bam" /dev/null
        samtools index -b -@ ~{cpus} "${filename}.aligned.sorted.bam"
    fi
    >>>

    output {
        File aligned_bam = basename(bam, ".bam") + ".aligned.sorted.bam"
        File aligned_bai = basename(bam, ".bam") + ".aligned.sorted.bam.bai"
        Int aligned_count = read_int("aligned_count.txt")
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
    filename=$(basename "~{bam}" .bam)
    samtools view --threads ~{cpus} -H "~{bam}" | grep "^@RG" > original_rg.txt

    samtools fastq --threads ~{cpus} -T "*" "~{bam}" | minimap2 -y --MD -ax map-ont -t ~{cpus} "~{ref_transcriptome}" - | \
    samtools sort --threads ~{cpus} | \
    samtools view --threads ~{cpus} -bh -F 260 -o temp.bam

    samtools view --threads ~{cpus} -H temp.bam > new_header.sam
    cat original_rg.txt >> new_header.sam
    samtools reheader new_header.sam temp.bam > "${filename}.transcriptome.aligned.sorted.bam"
    rm temp.bam new_header.sam original_rg.txt
    samtools index -b -@ ~{cpus} "${filename}.transcriptome.aligned.sorted.bam"
    if [ -s "${filename}.transcriptome.aligned.sorted.bam" ]; then
        echo "Alignment successful for ${filename}.bam"
        samtools view --threads ~{cpus} -c "${filename}.transcriptome.aligned.sorted.bam" > aligned_count.txt
        
    else
        echo "Alignment failed for ${filename}.bam, creating empty BAM file."
        echo 0 > aligned_count.txt
        samtools view -hb -o "${filename}.transcriptome.aligned.sorted.bam" /dev/null
        samtools index -b -@ ~{cpus} "${filename}.transcriptome.aligned.sorted.bam"
    fi
    >>>

    output {
        File aligned_bam = basename(bam, ".bam") + ".transcriptome.aligned.sorted.bam"
        File aligned_bai = basename(bam, ".bam") + ".transcriptome.aligned.sorted.bam.bai"
        Int aligned_count = read_int("aligned_count.txt")
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
    filename=$(basename "~{genome_bam}" .bam)
    NanoComp --bam "~{genome_bam}" "~{transcriptome_bam}" \
    --names "genome" "transcriptome" \
    --threads ~{cpus} \
    --outdir "nanocomp_report_${filename}"

    tar -czf "nanocomp_report_${filename}.tar.gz" "nanocomp_report_${filename}"
    >>>

    output {
        File report_tar = "nanocomp_report_" + basename(genome_bam, ".bam") + ".tar.gz"
    }

    runtime {
        cpu: cpus
        memory: "128GB"
        maxRunTime: 43200 #12 hours (12 * 3600 seconds)
        runtime_minutes: 1440 #24 hours (24 * 60 minutes)
        docker: "luxendr13/nanocomp:0.6.0"
        failOnStderr: false
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
    filename=$(basename "~{bam}" .bam)
    modkit pileup \
        ~{bam} \
        "${filename}.genome.bed" \
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
        --log-filepath "${filename}.genome.log" \
        --bedrmod
    >>>

    output {
        File bed = basename(bam, ".bam") + ".genome.bed"
        File log = basename(bam, ".bam") + ".genome.log"
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
    filename=$(basename "~{bam}" .bam)
    modkit pileup \
        ~{bam} \
        "${filename}.transcriptome.bed" \
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
        --log-filepath "${filename}.transcriptome.log" \
        --preload-references \
        --bedrmod
    >>>

    output {
        File bed = basename(bam, ".bam") + ".transcriptome.bed"
        File log = basename(bam, ".bam") + ".transcriptome.log"
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
