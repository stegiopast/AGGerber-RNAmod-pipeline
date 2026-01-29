## Usage

Run the pipeline by pointing the arguments to your local data and reference files.

**GPU & Model Configuration**
You can adjust the model parameters based on your available GPU resources and accuracy requirements:

* **Basecalling (`--model`):** Defaults to `sup` (Super High Accuracy). Switch to `hac` or `fast` to reduce GPU memory usage and increase processing speed.
* **Modifications (`--mods`):** Specifies which modifications to detect. You can detect all contexts (e.g., `m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG`) or a specific subset. Reduce the number of modifications to lower the computational load on your GPU.

```bash
nextflow run main.nf \
  --pod5_dir ./data/pod5_files/ \
  --model sup \
  --mods "m5C_2OmeC" \
  --ref_genome ./refs/genome.fa \
  --ref_transcriptome ./refs/transcripts.fa \
  --genome_gtf ./refs/annotation.gtf \
  --cpus 8
