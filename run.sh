nextflow run main.nf \
  --pod5_dir /home/johannes/100k_test_pod5/ \
  --model sup \
  --mods "m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG" \
  --ref_genome /home/johannes/RefGenomes/GRCh38.primary_assembly.genome.fa \
  --ref_transcriptome /home/johannes/RefGenomes/gencode.v47.transcripts.fa \
  --cpus 12
  
