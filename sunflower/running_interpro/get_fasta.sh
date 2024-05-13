#!/bin/bash
#SBATCH --job-name=get_fasta
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=40:00:00
#SBATCH --output=/scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/ha412_inflo_transcript_b3.fasta
#SBATCH --error=get_fasta.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL

#load samtools
ml SAMtools/1.10-GCC-8.3.0
# read in BEDTools
ml BEDTools/2.29.2-GCC-8.3.0


bedtools getfasta -fi /scratch/ely67071/sunflower_inflo_dev_data_b3/genomes/Ha412HO_v2.fasta -bed /scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/getfasta_bed_input_b3.bed -name
