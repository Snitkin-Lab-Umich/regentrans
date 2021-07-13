#subset Penn dataset for github
library(ape)

#load big data files ###############

#data is located in regentrans/extras/data_raw in a gzipped (.gz) format
#to use the raw data, first unzip the file(s) you will use using the command line
#to unzip a single file use the command: gunzip /regentrans/extras/data_raw/[desired file name].gz
#to unzip all files in the directory use the command: gunzup regentrans/extras/data_raw/*.gz

#metadata path Penn
metadata <- read.csv("/regentrans/extras/data_raw/ltach-metadata.csv")

#alignment path Penn
fasta <- ape::read.dna("/regentrans/extras/data_raw/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.filtered_polymorphic_sites.fasta",
                       format = "fasta")

#tree path Penn
tr <- ape::read.tree("/regentrans/extras/data_raw/2021_02_12_08_34_28_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

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
setwd("/regentrans/data")
saveRDS(Penn_test_input, file = "Penn_test_input_3.rds")

#test read back in
Penn_test_input2 <- readRDS(file = "/regentrans/data/Penn_test_input_3.rds")
locs2 <- Penn_test_input$locs
pt2 <- Penn_test_input$pt
fasta2 <- Penn_test_input$fasta
dists2 <- Penn_test_input$dists
tr2 <- Penn_test_input$tr

##save penn data
Penn_all_input <- list("locs" = locs_Penn, "pt" = pt_Penn, "fasta" = fasta_Penn, "dists" = dists_Penn, "tr" = tr_Penn)

getwd()
setwd("/regentrans/extras")
saveRDS(Penn_all_input, file = "Penn_test_input_all.rds")

Penn_all_input <- readRDS(file = "/regentrans/extras/Penn_test_input_all.rds")
locs_Penn <- Penn_all_input$locs
pt_Penn <- Penn_all_input$pt
fasta_Penn <- Penn_all_input$fasta
dists_Penn <- Penn_all_input$dists
tr_Penn <- Penn_all_input$tr

####save in RData file format
save(locs, file="/regentrans/data/locs.RData")
save(pt, file="/regentrans/data/pt.RData")
save(fasta, file="/regentrans/data/fasta.RData")
save(dists, file="/regentrans/data/dists.RData")
save(tr, file="/regentrans/data/tr.RData")


##################3/30/21 full data subset##################
library(dplyr)

#st info
st <- read.csv("regentrans/extras/data_raw/penn_kleborate_results.txt", sep = "\t")

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
save(locs, file="/regentrans/data/locs.RData")
save(pt, file="/regentrans/data/pt.RData")
save(dates, file="/regentrans/data/dates.RData")
save(fasta, file="/regentrans/data/fasta.RData")
save(dists, file="/regentrans/data/dists.RData")
save(tr, file="/regentrans/data/tr.RData")

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
save(pt_flow, file="/regentrans/data/pt_flow.RData")


####4/5/21
Fsp <- read.csv("Penn_Fsp_output.csv")
rownames(Fsp) <- Fsp[,1]
Fsp <- Fsp[,-1]
save(Fsp, file="/regentrans/data/Fsp.RData")
