#!/bin/bash
#SBATCH --job-name=convert_gff3_2_gtf
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=10GB
#SBATCH --time=10:00:00
#SBATCH --output=log_files/convertgff3.%j.out
#SBATCH --error=log_files/convert_gff3.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL

# november 2025


module load Cufflinks/20190706-GCC-12.3.0 

INPUT_DIR="$1"

gffread "$INPUT_DIR/Ha412HO_v2_braker.gff3" -T -o "${INPUT_DIR}/Ha412HO_v2_braker.gtf"
