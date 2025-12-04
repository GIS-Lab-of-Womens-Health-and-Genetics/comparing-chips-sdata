#!/usr/bin/env bash

plink='TODO/YOUR/PATH/TO/PLINK/plink.exe'
$plink --version

plink2='TODO/YOUR/PATH/TO/PLINK2/plink2.exe'
$plink2 --version

proj_dir=$(realpath "$(dirname $0)/..")
cd $proj_dir
echo "Running $0 from $(pwd)"

#### 1. Apply QC filters to GDA data
cd data/genotyping_arrays

$plink \
    --bfile data/genotyping_arrays/GDA_2025-02-17_clean \
    --geno 0.05 --maf 0.01 --hwe 0.0000001 --snps-only just-acgt \
    --make-bed \
    --out GDA_QC # for PCA 1000G

$plink \
    --bfile GDA_QC \
    --recode \
    --out GDA_QC_recode # for radmixture K7b dodecad


#### 2. Download 1000G phase 3 build 37 from $plink 2.0 resources - https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg
cd ../plink_pca_1000g

# rename psam
mv phase3_corrected.psam all_phase3.psam

# decompress pgen file
$plink2 --zst-decompress all_phase3.pgen.zst all_phase3.pgen

# convert to $plink 1 binary, removing multiallelic variants (not supported by bfiles)
$plink2 --pfile all_phase3 vzs --max-alleles 2 --make-bed --out 1000g_phase3_b37

# select subset of common SNPs from 1000G data, based on GDA SNPs after QC
$plink2 --bfile 1000g_phase3_b37 --extract ../genotyping_arrays/GDA_QC.bim --make-bed --geno 0.05 --maf 0.001 --chr 1-22 --out 1000g_GDA_subset


#### 3. Merge 1000G and GDA_QC data (subset of common biallelic variants)
# initial merge: creates .missnp file of variants with >2 alleles
$plink --bfile 1000g_GDA_subset --bmerge ../genotyping_arrays/GDA_QC --make-bed --out 1000GxGDA

# exclude variants with >2 alleles from 1000G and GDA data
$plink --bfile 1000g_GDA_subset --exclude 1000GxGDA-merge.missnp --make-bed --out temp-1000G
$plink --bfile ../genotyping_arrays/GDA_QC --exclude 1000GxGDA-merge.missnp --make-bed --out temp-GDA_QC

# remerge data after removing variants with 3+ alleles
$plink --bfile temp-1000G --bmerge temp-GDA_QC --make-bed --out 1000GxGDA


#### 4. Remove related samples (prune LD before IBD analysis)
# prune SNPs by LD: keep only "independent" SNPs
# pruned subset of markers that are in approximate linkage equilibrium with each other
# window size, step size, r^2 threshold
$plink2 --bfile 1000GxGDA --indep-pairwise 50 5 0.2 --out 1000GxGDA

# Identity-by-descent (IBD) on autosomal SNPs (after LD pruning) per pair of samples
$plink --bfile 1000GxGDA --genome gz --out 1000GxGDA --maf 0.01 --extract 1000GxGDA.prune.in

# create list of related samples
Rscript list_related_samples_1000G.R

# remove IBD-related samples
$plink --bfile 1000GxGDA --extract 1000GxGDA.prune.in --remove related_samples.txt --make-bed --out 1000GxGDA_final


#### 5. PCA
$plink2 --bfile 1000GxGDA_final --pca --out 1000GxGDA_pca

# plot PCA results
Rscript plot_PCA_results.R
