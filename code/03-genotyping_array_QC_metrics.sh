#!/usr/bin/env bash

plink='TODO/YOUR/PATH/TO/PLINK/plink.exe'
$plink --version

proj_dir=$(realpath "$(dirname $0)/..")
cd $proj_dir
echo "Running $0 from $(pwd)"

cd data/genotyping_arrays

# Get PLINK statistics
geno_arrs=( "GDA" "GSA" "OncoArray" "Axiom" )
for arr in "${geno_arrs[@]}"; do
  infile="${data_dir}/${arr}_${data_date}"
  out_prefix="$save_dir/${arr}_QC_$DATESTR"

  $plink --bfile ${arr}_2025-02-17_clean \
    --freq --missing --hardy --check-sex \
    --out ${arr}_stats
done

# Plot QC stats for all arrays
Rscript plot_QC_stats.R
