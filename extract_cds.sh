#!/bin/bash
#SBATCH --job-name=extract_CDS
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=15G
#SBATCH --time=8:00:00
#SBATCH --output=ErrorFiles/%x_%j.out
#SBATCH --error=ErrorFiles/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ely67071@uga.edu

ml Cufflinks/20190706-GCC-12.3.0

input_dir=$1




gffread "${input_dir}/lettuce_dev_data/genomes/Hildegardhap1v1_2.gff" -g "${input_dir}/lettuce_dev_data/genomes/Hildegardhap1v1.fasta" -x "lettuce/Hildegardhap1v1_cds.fasta"

gffread "${input_dir}/sunflower_dev_data/genomes/Ha412HO_v2_braker.gff3" -g "${input_dir}/sunflower_dev_data/genomes/Ha412HO_v2.fasta" -x "sunflower/Ha412HO_v2_cds.fasta"








