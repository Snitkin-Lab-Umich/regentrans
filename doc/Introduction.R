## ---- eval=FALSE, message = FALSE---------------------------------------------
#  # install devtools if you don't already have it installed
#  # install.packages('devtools')
#  
#  # install regentrans
#  devtools::install_github('Snitkin-Lab-Umich/regentrans')
#  
#  # load regentrans
#  library(regentrans)

## ----echo=FALSE, message=FALSE------------------------------------------------
devtools::load_all()

## ----setup,echo=TRUE, message=FALSE-------------------------------------------
library(ape)
library(tidyverse)
library(devtools)
library(ggtree)
library(pheatmap)
library(phytools)
library(gridExtra)
library(cowplot)

# set theme for plots 
theme_set(theme_bw() + theme(strip.background = element_rect(fill="white",linetype='blank'), text=element_text(size=15)))

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if your metadata is in a csv file
#  metadata <- readr::read_csv("/path/to/metadata.csv")

## ---- echo = FALSE, message = FALSE, include = FALSE--------------------------
pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if your alignment is in a fasta file
#  aln <- ape::read.dna("/path/to/aln.fasta",
#                       format = "fasta")

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if the tree is in Newick format
#  tr <- ape::read.tree("path/to/data.treefile")

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  dists <- ape::dist.dna(x = aln, # DNAbin object as read in above
#                         as.matrix = TRUE, # return as matrix
#                         model = "N", # count pairwise distances
#                         pairwise.deletion = TRUE # delete sites with missing data in a pairwise way
#                         )

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  dists <- ape::cophenetic.phylo(x = tr) # tree object as read in above

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if your patient transfer network is in a csv file
#  pt_trans_df <- readr::read_csv("/path/to/pt_trans_net.csv")

## -----------------------------------------------------------------------------
head(metadata)

## -----------------------------------------------------------------------------
class(aln)
aln

## -----------------------------------------------------------------------------
dists[1:5,1:5]

## -----------------------------------------------------------------------------
tr

## -----------------------------------------------------------------------------
head(pt_trans_df)

## -----------------------------------------------------------------------------
# named vector of locations
locs <- metadata %>% select(isolate_id, facility) %>% deframe()
head(locs)

## -----------------------------------------------------------------------------
ref_genome_length <- 5394056
mutation_rate <- 1.03e-6 # estimated from a K. pneumoniae ST258 time tree
floor(2*ref_genome_length*mutation_rate)

## -----------------------------------------------------------------------------
# get pair types for pairwise SNV distances (intra vs. inter)
pair_types <- get_pair_types(dists = dists, locs = locs, pt = pt)

# get fraction of intra-facility pairs for each SNV distance
frac_intra <- get_frac_intra(pair_types = pair_types)

## -----------------------------------------------------------------------------
head(pair_types)

## -----------------------------------------------------------------------------
head(frac_intra)

## -----------------------------------------------------------------------------
# plot fraction of intra-facility pairs for each SNV distance
frac_intra %>% 
  mutate(under = ifelse(pairwise_dist <= 10, '≤ 10 SNVs','Over threshold'), under = ifelse(pairwise_dist <= 6, '≤ 6 SNVs',under),
         under = factor(under, levels = c('≤ 6 SNVs', '≤ 10 SNVs','Over threshold'))) %>% 
  ggplot(aes(x = pairwise_dist, y = frac_intra, fill = under)) + 
  geom_bar(stat = "identity", alpha = 0.5) + 
  scale_fill_grey() + 
  labs(x = "Pairwise SNV distance", y = "Fraction of intra-facility pairs", fill = 'Possible SNV\nthresholds') + 
  ylim(0, 1) + xlim(-1,51) 

## -----------------------------------------------------------------------------
frac_intra_bin <- frac_intra %>% 
  filter(pairwise_dist <= 50) %>% 
  mutate(bin = cut(pairwise_dist, seq(0, max(pairwise_dist) + 5, 5), right = FALSE)) %>% 
  group_by(bin) %>% 
  summarise(n_intra = sum(n_intra), n_inter = sum(n_inter)) %>% 
  mutate(frac_intra = ifelse(n_intra != 0,
                             n_intra/(n_intra+n_inter), 0))


frac_intra_bin %>% 
  ggplot(aes(x = bin, y = frac_intra)) + 
  geom_bar(stat = "identity", alpha = 0.5) + 
  scale_fill_grey() + 
  labs(x = "Pairwise SNV distance", y = "Fraction of intra-facility pairs") + 
  ylim(0, 1) 

## ---- message=FALSE-----------------------------------------------------------
# pick prettier colors than the default
cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000')
names(cols) <- unique(locs)

# create vector for plotting location as colors
locs_tip <- c(locs[tr$tip.label], rep(NA, Nnode(tr)))

# plot tree
ggtree(tr) + geom_tippoint(aes(col = locs_tip)) + 
  scale_color_manual(values = cols) +
  labs(col = 'Facility') 

## ---- message=FALSE-----------------------------------------------------------
# get clusters (note this takes some time to run, but not too long for the test dataset!)
clusters <- get_clusters(tr,locs, pureness = 1, pt = pt)

## -----------------------------------------------------------------------------
# dissect output
pure_subtree_info <- clusters$pure_subtree_info
subtrees <- clusters$subtrees
pure_subtree_info

## -----------------------------------------------------------------------------
ggplot(data = pure_subtree_info, aes(x = loc, y = subtr_size, color = loc)) + 
  geom_jitter(position = position_jitter(width = 0.2, height = 0.1), alpha = 0.5) + 
  scale_color_manual(values = cols) +
  labs(y = "Number of isolates in \nphylogenetic cluster", x = "", color = 'Facility') +
  ylim(c(0,12.5)) +
  coord_flip()

## -----------------------------------------------------------------------------
# get ids of isolates in each cluster
isolates_clusters <- sapply(1:nrow(pure_subtree_info), function(x){
  i <- pure_subtree_info$index[x]
  name <- pure_subtree_info$isolate_id[x]
  if(!is.na(i)){
    name <- subtrees[[i]]$tip.label
  }
  name
})

# get unique cluster names
facil_clusters <- sapply(unique(locs), function(x){
  clusts <- isolates_clusters[pure_subtree_info$loc == x]
  names(clusts) <- paste0(x, 1:length(clusts))
  cluster_nums <- unlist(sapply(names(clusts), function(x){
  rep(x, length(clusts[[x]]))
  }))
  names(cluster_nums) <- unlist(clusts)
  cluster_nums
}) %>% unname() %>% unlist()

head(facil_clusters)

## -----------------------------------------------------------------------------
clust_tips <- c(facil_clusters[tr$tip.label], rep(NA, Nnode(tr)))
ggtree(tr) + geom_tippoint(aes(col = clust_tips)) + 
  theme(legend.position = 'none') 

## ---- message=FALSE-----------------------------------------------------------
snv_hist <- pair_types %>% 
  ggplot(aes(x = pairwise_dist, fill = pair_type)) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 30) +
  labs(x = "Pairwise SNV distance", y = "Count", fill = "") +
  geom_vline(xintercept = 11, col = 'darkgrey', size = 1) +
  theme(legend.justification = "top", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())
snv_hist

## ---- message=FALSE-----------------------------------------------------------
snv_hist + xlim(c(0,50))

## ---- warning = FALSE, message = FALSE----------------------------------------
left_join(pair_types, days_between_isolates) %>% 
  ggplot(aes(x = days_diff, y = pairwise_dist)) +
  geom_bin_2d() +
  geom_smooth(method = 'lm', color = 'white') +
  labs(x = '# days between samples', y = 'Pairwise SNV distance', fill = '# pairs')

## ---- message = FALSE---------------------------------------------------------
# create toy data
# each patient and each environment is a different location
# get patients 
pts <- unique(locs)[1:round(length(unique(locs))/2)]
# get environments
envs <- unique(locs)[(round(length(unique(locs))/2)+1):length(unique(locs))]
pair_types %>% mutate(src1 = ifelse(loc1 %in% pts, 'pt','env'),
                     src2 = ifelse(loc2 %in% pts, 'pt','env')) %>% 
  ggplot(aes(x = pairwise_dist, fill = paste(src1, src2))) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 30) +
  labs(x = "Pairwise SNV distance", y = "Count", fill = "") +
  geom_vline(xintercept = 11, col = 'darkgrey', size = 1)

## -----------------------------------------------------------------------------
fsp[1:4,1:4]

## ---- eval=TRUE, message=FALSE------------------------------------------------
# this takes a while to run, so we've pre-computed it and included it in the package
# fsp <- get_genetic_flow(fasta = aln, locs = locs, pt = pt)
pheatmap(fsp)

## ---- message=FALSE-----------------------------------------------------------
# get summary of pairwise SNV distance data
snv_summary <- summarize_pairs(pair_types,summary_fns = c('min'), threshs = c(6, 10)) 
snv_summary

## ---- message=FALSE-----------------------------------------------------------
# because each pair is only represented once, we must duplicate pairs with different locations
snv_summary_2 <- snv_summary
colnames(snv_summary_2) <- c("loc2", "loc1", "dist_min", "leq_6", "leq_10")
snv_summary_2 <- snv_summary_2 %>% relocate(loc1, .before = loc2) %>% filter(loc1 != loc2)
snv_summary_comb <- rbind(snv_summary, snv_summary_2)

# get facility summary
# facil_summary <- snv_summary_comb %>% 
#     mutate(pair_type = factor(ifelse(loc1 == loc2, 'Intra-facility pair', 'Inter-facility pair'), 
#                             levels = c('Intra-facility pair', 'Inter-facility pair'))) %>% 
#   group_by(loc1, pair_type, leq_6, leq_10) %>% summarize() %>% 
#   pivot_longer(cols = c(leq_6, leq_10)) %>% rename(thresh = name, n_pairs = value) %>% 
#   mutate(thresh = paste(gsub('leq_','≤ ', thresh), 'SNVs')) 

facil_summary <- snv_summary_comb %>% 
    mutate(pair_type = factor(ifelse(loc1 == loc2, 'Intra-facility pair', 'Inter-facility pair'), 
                            levels = c('Intra-facility pair', 'Inter-facility pair'))) %>% 
  group_by(loc1, pair_type) %>% summarize(leq_6 = sum(leq_6), leq_10 = sum(leq_10)) %>% 
  pivot_longer(cols = c(leq_6, leq_10)) %>% rename(thresh = name, n_pairs = value) %>% 
  mutate(thresh = paste(gsub('leq_','≤ ', thresh), 'SNVs')) 


# get locations to keep (at least 1 closely related pair)
keep_locs <- facil_summary %>% group_by(loc1) %>% summarize(sum_all = sum(n_pairs)) %>% filter(sum_all != 0) %>% select(loc1) %>% deframe()

# plot
facil_summary %>% filter(loc1 %in% keep_locs) %>% 
  tidyr::complete(loc1, pair_type, thresh) %>% 
  ggplot(aes(x = reorder(loc1, -n_pairs, sum), y = n_pairs, fill = pair_type)) +
  geom_col(position = 'dodge') +
       scale_fill_manual(values=c("cadetblue", "salmon")) + 
         labs(y = "Number of closely related \nisolate pairs", x = "Facility", fill = '', color = '') +
  facet_grid(factor(thresh, levels = c('≤ 6 SNVs','≤ 10 SNVs'))~.)

## ---- message=FALSE-----------------------------------------------------------
# get dates 
dates <- metadata %>% select(isolate_id, collection_date) %>% deframe()

time1 <- pair_types %>% 
  mutate(date1 = dates[isolate1], date2 = dates[isolate2]) %>% 
  filter(date1 == date2) %>% 
  ggplot(aes(x = pairwise_dist, fill = pair_type)) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 30) +
  labs(x = "Pairwise SNV distance", y = "Count", fill = "") +
  facet_grid(~date1) + 
  theme(legend.justification = "top", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

time2 <- pair_types %>% 
  mutate(date1 = dates[isolate1], date2 = dates[isolate2]) %>% 
  filter(date1 == date2) %>% 
  ggplot(aes(x = pairwise_dist, fill = pair_type)) + 
  geom_histogram(position = position_fill(), alpha = 0.4, bins = 30) +
  # xlim(0,20) +
  labs(x = "Pairwise SNV distance", y = "Fraction", fill = "") +
  facet_grid(~date1) + 
  theme(legend.justification = "top", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

plot_grid(time1, time2, nrow = 2, labels = 'AUTO')

## -----------------------------------------------------------------------------
head(pt_trans_df)

## -----------------------------------------------------------------------------
pt_trans <- get_patient_flow(pt_trans_df = pt_trans_df)
head(pt_trans)

## -----------------------------------------------------------------------------
# have to convert fsp matrix to long form first
fsp_long <- make_long_form(fsp)

# only keep inter-facility info in snv_summary
inter_snv_summary <- snv_summary %>% filter(loc1 != loc2)


pair_info <- merge_inter_summaries(pt_trans, inter_snv_summary, fsp_long) 
pair_info

## -----------------------------------------------------------------------------
pair_info %>% 
  ggplot(aes(x=sum_pt_trans_metric,y=fsp)) +
  geom_point() + geom_smooth(method='lm') + scale_x_log10() +
  labs(x = 'Patient flow', y = 'Genetic flow (Fsp)') 

## -----------------------------------------------------------------------------
pair_info %>% 
  ggplot(aes(x=sum_pt_trans_metric,y=leq_10)) +
  geom_point() + geom_smooth(method='lm') + scale_x_log10() +
  labs(x = 'Patient flow', y = '# closely related\npairs (≤ 10 SNVs)') 

## -----------------------------------------------------------------------------
pair_info %>%
  ggplot(aes(x=sum_pt_trans_metric,y=dist_min)) +
  geom_point() + geom_smooth(method='lm') + scale_x_log10() +
  labs(x = 'Patient flow', y = 'Minimum pairwise SNV distance') + ylim(c(0,50))

## -----------------------------------------------------------------------------
facil_coord

## -----------------------------------------------------------------------------
# get counts for each location
loc_n <- locs %>% table() %>% as.data.frame() %>% `colnames<-`(c("facil", "n_isolates"))

# merge longitude/latitude with locations
facil_geo <- facil_coord %>% left_join(loc_n, by = 'facil')
facil_geo

## -----------------------------------------------------------------------------
# merge facil_info and facil_geo
snv_geo <- pair_info %>% 
  left_join(facil_geo, by = c("loc1" = "facil")) %>% 
  rename(facil_1_lat = lat, facil_1_long = long, facil_1_n = n_isolates) %>% 
  left_join(facil_geo, by = c("loc2" = "facil")) %>% 
  rename(facil_2_lat = lat, facil_2_long = long, facil_2_n = n_isolates) %>% 
  filter(!is.na(facil_1_lat) & !is.na(facil_2_lat))

# filter to not plot zeros
snv_geo_thresh <- snv_geo %>% filter(leq_6 > 0) 

facil_geo %>% 
  ggplot(aes(x=long, y=lat)) +
  geom_curve(aes(x = facil_1_long, y = facil_1_lat, xend = facil_2_long, yend = facil_2_lat,
                 alpha = leq_6),
             data = snv_geo_thresh, curvature = 0.33) +
  geom_point(aes(size=n_isolates), alpha=1, color = "lightblue") + 
  theme_map() + 
  labs(size = 'Number of samples', alpha = 'Number of closely related\npairs (≤ 6 SNVs)')

## -----------------------------------------------------------------------------
snv_geo <- snv_geo %>% 
  mutate(geo_dist = sqrt((facil_1_long-facil_2_long)^2+(facil_1_lat-facil_2_lat)^2))

## -----------------------------------------------------------------------------
snv_geo %>% 
  ggplot(aes(x = geo_dist, y = leq_10)) + geom_point() + geom_smooth(method = 'lm') +
  labs(x = 'Physical distance between facilities', y = '# closely related pairs\n(≤ 10 SNVs)')

## -----------------------------------------------------------------------------
# merge pair_types and metadata, and find time difference
snv_and_metadat <- left_join(pair_types, metadata, by = c('isolate1' = 'isolate_id')) %>% 
  rename(patient_id_1 = patient_id, collection_date_1 = collection_date) %>% 
  left_join(metadata, by = c('isolate2' = 'isolate_id')) 
head(snv_and_metadat)

## ---- eval = FALSE------------------------------------------------------------
#  doFuture::registerDoFuture()
#  future::plan(future::multisession, workers = 4)
#  
#  clusters <- get_clusters(tr, locs)

