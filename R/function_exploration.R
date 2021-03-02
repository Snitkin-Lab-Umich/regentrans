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
#source get_largest_subtree
source("/Users/sophiehoffman/Desktop/regentrans/R/get_largest_subtree.R")
#source reverse_list_str
source("/Users/sophiehoffman/Desktop/regentrans/R/reverse_list_str.R")
#source get_clusters
source("/Users/sophiehoffman/Desktop/regentrans/R/get_clusters.R")
#source get_facility_fsp
source("/Users/sophiehoffman/Desktop/regentrans/R/get_facility_fsp.R")


#devtools::load_all("/Users/sophiehoffman/Desktop/regentrans")
#load_all()

###################realm data prep#######################
#metadata path Realm
#/nfs/turbo/umms-esnitkin/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv
#metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv")

#alignment path Realm
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
#fasta <- ape::read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta",
#                       format = "fasta")

#tree path Realm
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile
#tree <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile")

#locs - locations of isolates (e.g. facility of isolation)
#named vector, where names are same as names of dist.dna output
#locs <- metadata$f_id
#names(locs) <- metadata$gID

#pt
#pt <- metadata$pt_id
#names(pt) <- metadata$gID

#for the test data I have to remove the last char of each of the row and col names
#rownames(dists) <- substr(rownames(dists), 1, nchar(rownames(dists))-1)
#colnames(dists) <- substr(colnames(dists), 1, nchar(colnames(dists))-1)

##########################################################
###################penn data prep#########################
#metadata path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv
#metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv")

#alignment path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
#fasta <- ape::read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta",
                       format = "fasta")

#tree path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile
#tree <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

#locs - locations of isolates (e.g. facility of isolation)
#named vector, where names are same as names of dist.dna output
#locs <- metadata$ltach
#names(locs) <- paste0("PCMP_H", metadata$isolate_no)

#pt
#pt <- metadata$patient_id
#names(pt) <- paste0("PCMP_H", metadata$isolate_no)
##########################################################
##################github data prep
##########################################################
Penn_test_input <- readRDS(file = "/Users/sophiehoffman/Desktop/regentrans/extras/Penn_test_input.rds")
locs <- Penn_test_input$locs
pt <- Penn_test_input$pt
fasta <- Penn_test_input$fasta
dists <- Penn_test_input$dists

#############################################################################################################################################
#prep for the get_snv_dists function

#dists - snv distance matrix returned by dist.dna
#use as.matrix = true for the dist.dna function
#dists <- dist.dna(x = fasta, as.matrix = TRUE, model = "N")

#run the function
snv_dists <- get_snv_dists(dists, locs, pt)
snv_dist_no_pt <- get_snv_dists(dists, locs)

#for testing purposes, subset to ones that are CRE only
#maybe we will use this in the uploaded thing? or will we be using penn data? if so what does that look like?
#isolates <- intersect(names(locs), rownames(dists))
#subset the DNAbin object to the samples they have in common
#dists<-dists[isolates,isolates]
#############################################################################################################################################
#get_frac_intra
threshs <- seq(1,19041, by=1)
frac_intra <- get_frac_intra(snv_dists = snv_dists, threshs = threshs)
frac_intra_2 <- get_frac_intra(dists = dists, locs = locs, pt = pt, threshs = threshs)

#############################################################################################################################################
#get_clusters
clusters <- get_clusters(tr,locs)
#dissect output
pure_subtree_info <- clusters$pure_subtree_info
subtrees <- clusters$subtrees
cluster_pureness <- clusters$cluster_pureness


#############################################################################################################################################
#get_facility_fsp
fsp <- get_facility_fsp(fasta,locs)

