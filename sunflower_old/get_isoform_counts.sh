#!/bin/bash
#!/bin/bash
#SBATCH --job-name=get_isoform_counts
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=100gb
#SBATCH --time=40:00:00


# written 10-17-2024


# Input parameters
gtf_directory="$1" # Directory containing the GTF files
output_prefix="$2" # Output prefix for gene and transcript counts

# Check if required arguments are provided
if [[ -z "$gtf_directory" || -z "$output_prefix" ]]; then
  echo "Usage: $0 <gtf_directory> <output_prefix>"
  exit 1
fi

# Create output files
transcript_output="${output_prefix}_transcript_counts.csv"
gene_output="${output_prefix}_gene_counts.csv"

# Initialize CSV headers
echo -e "Transcript_ID,Gene_ID,Sample,FPKM,TPM,Cov" > "$transcript_output"
echo -e "Gene_ID,Sample,FPKM,TPM,Cov" > "$gene_output"

# Process each GTF file in the directory
for gtf_file in "$gtf_directory"/*.gtf; do
  sample=$(basename "$gtf_file" .gtf)  # Extract sample name from filename

  echo "Processing $sample..."

  # Extract transcript-level counts (FPKM, TPM, Cov) from each GTF file
  awk -v sample="$sample" '
    $3 == "transcript" {
      for(i=1; i<=NF; i++) {
        if($i ~ /transcript_id/) transcript_id=$(i+1)
        if($i ~ /gene_id/) gene_id=$(i+1)
        if($i ~ /FPKM/) fpkm=$(i+1)
        if($i ~ /TPM/) tpm=$(i+1)
        if($i ~ /cov/) cov=$(i+1)
      }
      gsub(/"|;/, "", transcript_id)
      gsub(/"|;/, "", gene_id)
      gsub(/"|;/, "", fpkm)
      gsub(/"|;/, "", tpm)
      gsub(/"|;/, "", cov)
      print transcript_id "," gene_id "," sample "," fpkm "," tpm "," cov
    }
  ' "$gtf_file" >> "$transcript_output"

  # Extract gene-level counts by summing transcript FPKM, TPM, and Cov for each gene
  awk -v sample="$sample" '
    $3 == "transcript" {
      for(i=1; i<=NF; i++) {
        if($i ~ /gene_id/) gene_id=$(i+1)
        if($i ~ /FPKM/) fpkm=$(i+1)
        if($i ~ /TPM/) tpm=$(i+1)
        if($i ~ /cov/) cov=$(i+1)
      }
      gsub(/"|;/, "", gene_id)
      gsub(/"|;/, "", fpkm)
      gsub(/"|;/, "", tpm)
      gsub(/"|;/, "", cov)
      gene_fpkm[gene_id] += fpkm
      gene_tpm[gene_id] += tpm
      gene_cov[gene_id] += cov
    }
    END {
      for(gene_id in gene_fpkm) {
        print gene_id "," sample "," gene_fpkm[gene_id] "," gene_tpm[gene_id] "," gene_cov[gene_id]
      }
    }
  ' "$gtf_file" >> "$gene_output"

done



