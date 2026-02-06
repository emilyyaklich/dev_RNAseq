#!/bin/bash
#SBATCH --job-name=check_strand
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=15G
#SBATCH --time=8:00:00
#SBATCH --output=log_files/%x_%j.out
#SBATCH --error=log_files/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ely67071@uga.edu

# check strandedness of RNAseq data
# 02-05-2026
ml RSeQC/5.0.4-foss-2023a 
ml BEDOPS/2.4.41-foss-2023a


input_dir=$1

gtf2bed < "${input_dir}/genomes/Hildegardhap1v1_2.gtf" > "${input_dir}/genomes/Hildegardhap1v1_2.bed"



infer_experiment.py -i "${input_dir}/STAR_output_gtf/CM4Aligned.sortedByCoord.out.bam" -r "${input_dir}/genomes/Hildegardhap1v1_2.bed"


