#walkthrough of functions, data loading

#load libraries
library(ape)
library(dplyr)
library(devtools)

#source checks
source("/Users/sophiehoffman/Desktop/regentrans/R/checks.R")
#source tests
source("/Users/sophiehoffman/Desktop/regentrans/R/tests.R")
#source get_snv_dists
source("/Users/sophiehoffman/Desktop/regentrans/R/get_snv_dists.R")
#source get_frac_intra
source("/Users/sophiehoffman/Desktop/regentrans/R/get_frac_intra.R")
#source get_frac_intra_from_snv_dists
source("/Users/sophiehoffman/Desktop/regentrans/R/get_frac_intra_from_snv_dists.R")
#source get_largest_subtree
source("/Users/sophiehoffman/Desktop/regentrans/R/get_largest_subtree.R")
#source reverse_list_str
source("/Users/sophiehoffman/Desktop/regentrans/R/reverse_list_str.R")


#devtools::load_all("/Users/sophiehoffman/Desktop/regentrans")
#load_all()

#metadata path
#/nfs/turbo/umms-esnitkin/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv
metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv")

#alignment path
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
fasta <- ape::read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta", format = "fasta")


#tree path
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile
tree <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile")

#############################################################################################################################################
#prep for the get_snv_dists function

#dists - snv distance matrix returned by dist.dna
#use as.matrix = true for the dist.dna function
dists <- dist.dna(x = fasta, as.matrix = TRUE)

#locs - locations of isolates (e.g. facility of isolation)
#named vector, where names are same as names of dist.dna output
locs <- metadata$f_id
names(locs) <- metadata$gID

#pt
pt <- metadata$pt_id
names(pt) <- metadata$gID

#for the test data I have to remove the last char of each of the row and col names
rownames(dists) <- substr(rownames(dists), 1, nchar(rownames(dists))-1)
colnames(dists) <- substr(colnames(dists), 1, nchar(colnames(dists))-1)
#run the function
snv_dists <- get_snv_dists(dists, locs, pt)

#############################################################################################################################################
#get_frac_intra
threshs <- seq(0.001,0.0264350436487801,0.005)
frac_intra <- get_frac_intra(snv_dists, threshs)

#############################################################################################################################################
#wrapper get_frac_intra_from_snv_dists
#maybe this will be the only outward facing function
frac_intra_2 <- get_frac_intra_from_snv_dists(dists, locs, pt, threshs)

#############################################################################################################################################
#get_clusters
tr <- tree
tr$tip.label <- substr(tr$tip.label, 1, nchar(tr$tip.label)-1)

#for running tests
locs_sub <-locs_sub[!is.na(locs_sub)]

#############################################################################################################################################
#get_facility_fsp
snp_dist <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_Nursing_Home/Sequence_data/final_genomics_files/2019-07-19_Filtered_Pathways_VREfm", sep = " ", header = FALSE, skip = 1)
#what are the dimensions? Is that what the first line of the file is?
dim(snp_dist)

