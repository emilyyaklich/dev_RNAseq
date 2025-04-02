#!/bin/bash
#SBATCH --job-name=run_mapq
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100GB
#SBATCH --time=40:00:00
#SBATCH --output=mapq.%j.out
#SBATCH --error=mapq.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL




module load MultiQC/1.14-foss-2022a

multiqc /scratch/ely67071/sunflower_dev_data/STAR_output/all_log_files/ -o /scratch/ely67071/sunflower_dev_data/STAR_output/all_log_files/
