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
tr <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

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
#subset the tree using keep.tip
tr_subset <- keep.tip(tr,names(locs))

dists <- dist.dna(x = fasta_sub, as.matrix = TRUE, model = "N")

#make list
#Penn_test_input <- list("locs" = locs_sub, "pt" = pt_sub, "fasta" = fasta_sub, "dists" = dists)
Penn_test_input <- list("locs" = locs, "pt" = pt, "fasta" = fasta, "dists" = dists, "tr" = tr)

getwd()
setwd("/Users/sophiehoffman/Desktop/regentrans/data")
saveRDS(Penn_test_input, file = "Penn_test_input_3.rds")

#test read back in
Penn_test_input2 <- readRDS(file = "/Users/sophiehoffman/Desktop/regentrans/data/Penn_test_input_3.rds")
locs2 <- Penn_test_input$locs
pt2 <- Penn_test_input$pt
fasta2 <- Penn_test_input$fasta
dists2 <- Penn_test_input$dists
tr2 <- Penn_test_input$tr

##save penn data
Penn_all_input <- list("locs" = locs_Penn, "pt" = pt_Penn, "fasta" = fasta_Penn, "dists" = dists_Penn, "tr" = tr_Penn)

getwd()
setwd("/Users/sophiehoffman/Desktop/regentrans/extras")
saveRDS(Penn_all_input, file = "Penn_test_input_all.rds")

Penn_all_input <- readRDS(file = "/Users/sophiehoffman/Desktop/regentrans/extras/Penn_test_input_all.rds")
locs_Penn <- Penn_all_input$locs
pt_Penn <- Penn_all_input$pt
fasta_Penn <- Penn_all_input$fasta
dists_Penn <- Penn_all_input$dists
tr_Penn <- Penn_all_input$tr

####save in RData file format
save(locs, file="/Users/sophiehoffman/Desktop/regentrans/data/locs.RData")
save(pt, file="/Users/sophiehoffman/Desktop/regentrans/data/pt.RData")
save(fasta, file="/Users/sophiehoffman/Desktop/regentrans/data/fasta.RData")
save(dists, file="/Users/sophiehoffman/Desktop/regentrans/data/dists.RData")
save(tr, file="/Users/sophiehoffman/Desktop/regentrans/data/tr.RData")


##################3/30/21 full data subset##################
library(dplyr)

#metadata path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv
metadata <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Analysis/regentrans_data/2021-02-16_subset-data/data/ltach-metadata.csv")

#alignment path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta
fasta <- ape::read.dna("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta",
                       format = "fasta")

#tree path Penn
#/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile
tr <- ape::read.tree("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/2021_02_10_Penn_All_variant_calling/2021_02_12_08_34_28_core_results/gubbins/iqtree_masked_wga/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

#st info
st <- read.csv("/Users/sophiehoffman/Desktop/gl_mount/Project_Penn_KPC/Sequence_data/Reports/kleborate/penn_kleborate_results.txt", sep = "\t")

#join st and metadata
metadata$isolate_no <- paste0("PCMP_H", metadata$isolate_no)
metadata <- metadata %>% left_join(st, by = c("isolate_no" = "strain")) %>% filter(ST == "ST258") %>% select(isolate_no, patient_id, ltach, cx_date)
#subset date to year
metadata$cx_date <- substr(metadata$cx_date, 1, 4)

#make locs
locs <- metadata$ltach
names(locs) <- metadata$isolate_no
#make pt
pt <- metadata$patient_id
names(pt) <- metadata$isolate_no
#make date
dates <- as.factor(metadata$cx_date)
names(dates) <- metadata$isolate_no

#find common
common <- intersect(intersect(names(locs), rownames(fasta)), tr$tip.label)
locs_sub <- locs[names(locs) %in% common]
pt_sub <- pt[names(pt) %in% common]
dates_sub <- dates[names(dates) %in% common]
#subset fasta
fasta_sub<-fasta[common,]
#subset the tree using keep.tip
tr_sub <- keep.tip(tr,common)
#make dists
dists <- dist.dna(x = fasta_sub, as.matrix = TRUE, model = "N", pairwise.deletion = TRUE)

#rename them all
locs <- locs_sub
pt <- pt_sub
dates <- dates_sub
fasta <- fasta_sub
tr <- tr_sub

#save all as .RData
save(locs, file="/Users/sophiehoffman/Desktop/regentrans/data/locs.RData")
save(pt, file="/Users/sophiehoffman/Desktop/regentrans/data/pt.RData")
save(dates, file="/Users/sophiehoffman/Desktop/regentrans/data/dates.RData")
save(fasta, file="/Users/sophiehoffman/Desktop/regentrans/data/fasta.RData")
save(dists, file="/Users/sophiehoffman/Desktop/regentrans/data/dists.RData")
save(tr, file="/Users/sophiehoffman/Desktop/regentrans/data/tr.RData")

#make a pt flow function
library(Matrix)
x<-Matrix(sample(20:100, length(unique(locs))^2, replace=TRUE),length(unique(locs)))
diag(x) <- 0
rownames(x) <- unique(locs)
colnames(x) <- unique(locs)
pat_flow <- na.omit(data.frame(as.table(as.matrix(x))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)
pt_flow <- pat_flow
save(pt_flow, file="/Users/sophiehoffman/Desktop/regentrans/data/pt_flow.RData")


