#walkthrough of functions, data loading

#load libraries
library(ape)

#metadata path
#/nfs/turbo/umms-esnitkin/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv
metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Analysis/NDM_transmission/2020-12-16_manuscript-figures/data/kp_st147_metadata.csv")

#alignment path
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
fasta <- read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta", format = "fasta")


#tree path
#/nfs/turbo/umms-esnitkin/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile
tree <- read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_REALM/Sequence_data/output_files/2020_12_09_variant_calling_ST147/2020_12_11_12_49_38_core_results/gubbins/iqtree_masked_wga/2020_12_11_12_49_38_Kp46596_genome_aln_w_alt_allele_unmapped.treefile")

#prep for the get_snv_dists function

#dists - snv distance matrix returned by dist.dna
#use as.matrix = true for the dist.dna function
dists <- dist.dna(x = fasta, as.matrix = TRUE)

#locs - locations of isolates (e.g. facility of isolation)
#using pt_id here
#will we do pre-processing outside or inside of the function definition
#named list, where names are same as names of dist.dna output
locs <- metadata$f_id
names(locs) <- metadata$name
#unname(locs[hi])
#hi <- c("REALM_CRE_165", "REALM_CRE_164")

#pt
pt <- metadata$pt_id
names(pt) <- metadata$name
#unname(pt[hi])



