#!/bin/bash
#SBATCH --job-name=run_stringtie
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=100gb
#SBATCH --time=40:00:00
#SBATCH --output=../ErrorFiles/run_stringtie.%j.out
#SBATCH --error=../ErrorFiles/run_stringtie.%j.error

# written on 10-17-2024

ml StringTie

bam_dir="/scratch/ely67071/sunflower_inflo_dev_data_b3/STAR_output/BAM_files/extracted"

# List of BAM files
bam_files=(
    "$bam_dir/extracted_region_10D_REP1_ATTACTCG-TATAGCCT_S97_L007Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_20D_REP2_TCCGGAGA-ATAGAGGC_S98_L007Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_30D_REP2_CGCTCATT-CCTATCCT_S99_L007Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_35D_REP1_GAGATTCC-GGCTCTGA_S100_L007Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_10D_2_ACCTTGGC-GGCCTCAT_S39_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_10D_3_ATATCTCG-ATCTTAGT_S40_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_20D_2_GCGCTCTA-GCTCCGAC_S41_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_20D_3_AACAGGTT-ATACCAAG_S42_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_30D_2_GGTGAACC-GCGTTGGA_S43_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_30D_3_CAACAATG-CTTCACGG_S44_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_35D_2_TGGTGGCA-TCCTGTAA_S45_L002Aligned.sortedByCoord.out.bam"
    "$bam_dir/extracted_region_HA_35D_3_AGGCAGAG-AGAATGCC_S46_L002Aligned.sortedByCoord.out.bam"
)



# Loop through the BAM files and run stringtie
for bam_file in "${bam_files[@]}"; do
    # Extract the base name without the .bam extension
    base_name=$(basename "$bam_file" .bam)
    # Run stringtie
    # -G use braker3 annotation to guide isoforms
    stringtie "$bam_file" -G /scratch/ely67071/sunflower_inflo_dev_data_b3/genomes/Ha412HO_v2_braker.gff3  -o "$bam_dir/${base_name}_output.gtf"
done

# merge files to consolidate isoforms
stringtie --merge -o $bam_dir/merged_output.gtf $bam_dir/*_output.gtf


# quantify merged gtf
for bam_file in "${bam_files[@]}"; do
    base_name=$(basename "$bam_file" .bam)
    stringtie "$bam_file" -e -B -o "$bam_dir/${base_name}_quantified.gtf" -G "$bam_dir/merged_output.gtf"
done


