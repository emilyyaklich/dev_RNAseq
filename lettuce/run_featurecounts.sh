#!/bin/bash
#SBATCH --job-name=run_featurecounts
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100GB
#SBATCH --time=48:00:00
#SBATCH --output=log_files/run_featurecounts.%j.out
#SBATCH --error=log_files/run_featurecounts.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL




module load Subread/2.0.6-GCC-12.3.0



ANNOTATION_DIR="$1"
BAM_FILE_DIR="$2"
OUTPUT_DIR="$3"

featureCounts -T 8 -p -B -C -s 2 --tmpDir /scratch/ely67071/feature_tmp  -a "$ANNOTATION_DIR/Hildegardhap1v1_2.gtf" -o "$OUTPUT_DIR/lettuce_dev_fc_gtf.txt" "$BAM_FILE_DIR"/*.bam

# featureCounts parameters explanation:
# -p   → treat input as paired-end
# -B   → only count fragments where both mates are successfully aligned
# -C   → only count properly paired fragments (mates mapped to the same chromosome with correct orientation and distance)
# -T 8 → number of threads to use (8 in this case)
# -s 2 → stranded library (reverse stranded reads) from running infer_experiment.py got these results

