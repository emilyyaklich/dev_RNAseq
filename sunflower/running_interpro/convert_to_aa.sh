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

ml EMBOSS/6.6.0-foss-2021b
transeq -sformat pearson -frame 6 -sequence /scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/ha412_inflo_transcript_b3.fasta -outseq /scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/ha412_inflo_aa_b3.fasta -clean 
