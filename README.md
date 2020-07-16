# PerHapS: Paired End Reads based HAPlotyping for Sequencing
PerHapS is a new and simple approach to directly call haplotypes from short-read, paired-end WES data
Author: Jie Huang, MD, PhD, Department of Global Health, Peking University School of Public Health


The technical bottleneck in direct haplotype calling from short-read sequencing lies in the length of sequenced DNA fragments, often too short to include two or multiple variable nucleotide positions that define the haplotype of interest. Indeed, while sequencing reads length in UKBB WES data is 76bp, the two APOE SNPs (rs7412 and rs429358), defining the common APOE polymorphism, are located 138 bp apart. We pieced short reads by utilizing their labels to generate a composite haplotype longer than 138bp.

Steps:

#1. download UKB pre-phased genetic data

#2. download UKB WES data through a look, extract a certain genetic region such as APOE to save storage space

#3. run Perhaps.sh to piece together haplotypes from WES

#4. run Perhaps.R and more analyses to explore the haplotypes

