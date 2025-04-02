#!/bin/bash
#SBATCH --job-name=extract_gene_region
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=5gb
#SBATCH --time=40:00:00
#SBATCH --output=../ErrorFiles/extract_gene_region.%j.out
#SBATCH --error=../ErrorFiles/extract_gene_region.%j.error

# written 10-17-2024

module load SAMtools
# Directory containing the BAM files
bam_dir="/scratch/ely67071/sunflower_inflo_dev_data_b3/STAR_output/BAM_files"

# List of BAM files
files=("$bam_dir/10D_REP1_ATTACTCG-TATAGCCT_S97_L007Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_30D_2_GGTGAACC-GCGTTGGA_S43_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_30D_3_CAACAATG-CTTCACGG_S44_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/20D_REP2_TCCGGAGA-ATAGAGGC_S98_L007Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_35D_2_TGGTGGCA-TCCTGTAA_S45_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/30D_REP2_CGCTCATT-CCTATCCT_S99_L007Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_35D_3_AGGCAGAG-AGAATGCC_S46_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/35D_REP1_GAGATTCC-GGCTCTGA_S100_L007Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_10D_3_ATATCTCG-ATCTTAGT_S40_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_20D_2_GCGCTCTA-GCTCCGAC_S41_L002Aligned.sortedByCoord.out.bam" \
       "$bam_dir/HA_20D_3_AACAGGTT-ATACCAAG_S42_L002Aligned.sortedByCoord.out.bam"
      )


# Loop through each BAM file
for file in "${files[@]}"; do
  # Extract the file ID (remove the .bam extension)
  fileID=$(basename "$file" .bam)
  
  # Index the BAM file (necessary to extract)
  samtools index "$file"
  
  # Extract SPLAYED gene region
  samtools view -h "$file" "Ha412HOChr17:100804813-100824458" > "$bam_dir/extracted_region_${fileID}.sam"
  
  # Convert the SAM to BAM format (necessary for indexing)
  samtools view -bS "$bam_dir/extracted_region_${fileID}.sam" > "$bam_dir/extracted_region_${fileID}.bam"
  
  # BAM file
  samtools index "$bam_dir/extracted_region_${fileID}.bam"
  
done


