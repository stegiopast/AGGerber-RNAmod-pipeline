#!/usr/bin/env bash
set -euo pipefail

POD5_DIR="/home/johannes/100k_test_pod5/split_pod5"
SAMPLE_ID="ONT_mRNA_pilot_sample_TEST"
REF_GENOME="/home/johannes/RefGenomes/v49/GRCh38.primary_assembly.genome.fa"
REF_TX="/home/johannes/RefGenomes/v49/gencode.v49.transcripts.fa"
CPUS=12

POD5_JSON=$(printf '%s\n' "$POD5_DIR"/*.pod5 | jq -R . | jq -s .)

cat > inputs.json <<EOF
{
  "ont_mRNA_pilot.pod5_files": ${POD5_JSON},
  "ont_mRNA_pilot.sample_id": "${SAMPLE_ID}",
  "ont_mRNA_pilot.ref_genome": "${REF_GENOME}",
  "ont_mRNA_pilot.ref_transcriptome": "${REF_TX}",
  "ont_mRNA_pilot.cpus": ${CPUS}
}
EOF