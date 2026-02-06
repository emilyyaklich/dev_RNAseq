#!/bin/bash
#SBATCH --job-name=convert_gff_2_gtf
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=10GB
#SBATCH --time=10:00:00
#SBATCH --output=log_files/convertgff3.%j.out
#SBATCH --error=log_files/convert_gff3.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL

# November 2025


module load AGAT/1.4.2-GCC-13.3.0

INPUT_DIR="$1"

agat_convert_sp_gff2gtf.pl --gff "$INPUT_DIR/Hildegardhap1v1_2.gff" -o "${INPUT_DIR}/Hildegardhap1v1_2_v2.gtf"
