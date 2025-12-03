#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=15G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ely67071@uga.edu

ml Cufflinks/20190706-GCC-12.3.0

input_dir=$1


gffread "${input_dir}/genomes/Hildegardhap1v1_2.gff" -g "${input_dir}/genomes/Hildegardhap1v1.fasta" -y "${input_dir}/entap/lettuce_aa_Hildegardhap1v1.fasta"
