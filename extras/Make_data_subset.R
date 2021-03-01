#subset Penn dataset for github 
library(ape)

#load big data files ###############

#metadata path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv
metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv")

#alignment path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
fasta <- ape::read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta",
                       format = "fasta")

#tree path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile
tree <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

#locs - locations of isolates (e.g. facility of isolation)
#named vector, where names are same as names of dist.dna output
locs <- metadata$ltach
names(locs) <- paste0("PCMP_H", metadata$isolate_no)

#pt
pt <- metadata$patient_id
names(pt) <- paste0("PCMP_H", metadata$isolate_no)

#subset big data files 
locs_sub <- locs[1:102]
pt_sub <- pt[1:102]
common <- intersect(names(locs_sub), rownames(fasta))
locs_sub <- locs[names(locs) %in% common]
pt_sub <- pt[names(pt) %in% common]
fasta_sub<-fasta[common,]

dists <- dist.dna(x = fasta_sub, as.matrix = TRUE, model = "N")

#make list 
Penn_test_input <- list("locs" = locs_sub, "pt" = pt_sub, "fasta" = fasta_sub, "dists" = dists)

getwd()
setwd("/Users/sophiehoffman/Desktop/regentrans/extras")
saveRDS(Penn_test_input, file = "Penn_test_input.rds")

#test read back in 
Penn_test_input <- readRDS(file = "/Users/sophiehoffman/Desktop/regentrans/extras/Penn_test_input.rds")
locs <- Penn_test_input$locs
pt <- Penn_test_input$pt
fasta <- Penn_test_input$fasta
dists <- Penn_test_input$dists
