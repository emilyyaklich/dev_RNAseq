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


gffread "${input_dir}/genomes/Ha412HO_v2_braker.gff3" -g "${input_dir}/genomes/Ha412HO_v2.fasta" -y Ha412HO_v2_aa_braker3.fasta
