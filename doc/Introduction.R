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

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if your alignment is in a fasta file
#  aln <- ape::read.dna("/path/to/aln.fasta",
#                       format = "fasta")

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  dists <- ape::dist.dna(x = aln, # DNAbin object as read in above
#                         as.matrix = TRUE, # return as matrix
#                         model = "N", # count pairwise distances
#                         pairwise.deletion = TRUE # delete sites with missing data in a pairwise way
#                         )

## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # this is if the tree is in Newick format
#  tr <- ape::read.tree("path/to/data.treefile")

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
# subset metadata to include the first isolate from each patient per facility
metadata <- metadata %>% arrange(collection_date) %>% filter(!duplicated(patient_id, facility))

# isolate subset
isolate_subset <- metadata$sample_id

# subset alignment
aln <- aln[isolate_subset,]

# subset tree
tr <- keep.tip(tr, isolate_subset)

# subset distance matrix
dists <- dists[isolate_subset, isolate_subset]

## -----------------------------------------------------------------------------
# named vector of locations
locs <- metadata %>% select(sample_id, facility) %>% deframe()
head(locs)

## -----------------------------------------------------------------------------
ref_genome_length <- 5394056
mutation_rate <- 1.03e-6 # estimated from a K. pneumoniae ST258 time tree
round(2*ref_genome_length*mutation_rate)

## -----------------------------------------------------------------------------
# get pair types for pairwise SNV distances (intra vs. inter)
snv_dists <- get_snv_dists(dists = dists, locs = locs)

# get fraction of intra-facility pairs for each SNV distance
frac_intra <- get_frac_intra(snv_dists = snv_dists)

## -----------------------------------------------------------------------------
head(snv_dists)

## -----------------------------------------------------------------------------
head(frac_intra)

## -----------------------------------------------------------------------------
# plot fraction of intra-facility pairs for each SNV distance
frac_intra %>% 
  mutate(under = ifelse(pairwise_dist < 11, '< 11 SNVs','Over threshold'), under = ifelse(pairwise_dist < 7, '< 7 SNVs',under),
         under = factor(under, levels = c('< 7 SNVs', '< 11 SNVs','Over threshold'))) %>% 
  ggplot(aes(x = pairwise_dist, y = frac_intra, fill = under)) + 
  geom_bar(stat = "identity", alpha = 0.5) + 
  scale_fill_grey() + 
  labs(x = "Pairwise SNV distance", y = "Fraction of intra-facility pairs", fill = 'Possible SNV\nthresholds') + 
  ylim(0, 1) + xlim(-1,51) 

## -----------------------------------------------------------------------------
frac_intra_bin <- frac_intra %>% 
  filter(pairwise_dist < 50) %>% 
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
clusters <- get_clusters(tr,locs, pureness = 1)

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
  name <- pure_subtree_info$sample_id[x]
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
snv_hist <- snv_dists %>% 
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

## ---- message = FALSE---------------------------------------------------------
# create toy data
# each patient and each environment is a different location
# get patients 
pts <- unique(locs)[1:round(length(unique(locs))/2)]
# get environments
envs <- unique(locs)[(round(length(unique(locs))/2)+1):length(unique(locs))]
snv_dists %>% mutate(src1 = ifelse(loc1 %in% pts, 'pt','env'),
                     src2 = ifelse(loc2 %in% pts, 'pt','env')) %>% 
  ggplot(aes(x = pairwise_dist, fill = paste(src1, src2))) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 30) +
  labs(x = "Pairwise SNV distance", y = "Count", fill = "") +
  geom_vline(xintercept = 11, col = 'darkgrey', size = 1)

## -----------------------------------------------------------------------------
fsp[1:4,1:4]

## ---- eval=FALSE, message=FALSE-----------------------------------------------
#  # this takes a while to run, so we've pre-computed it and included it in the package
#  # fsp <- get_facility_fsp(aln, locs)
#  pheatmap(fsp)

## ---- message=FALSE-----------------------------------------------------------
close_pair_dat <- bind_rows(bind_cols(snvthresh = '< 7 SNVs', snv_dists %>% filter(pairwise_dist < 7) %>%
                     mutate(loc1=factor(loc1), pair_type=factor(pair_type)) %>%
                     group_by(loc1, pair_type) %>% summarize(n = n())
                     %>% tidyr::complete(loc1, pair_type)),
                    bind_cols(snvthresh = '< 11 SNVs', snv_dists %>% filter(pairwise_dist < 11) %>%
                     mutate(loc1=factor(loc1), pair_type=factor(pair_type)) %>%
                     group_by(loc1, pair_type) %>% summarize(n = n())
                     %>% tidyr::complete(loc1, pair_type)))

close_pair_dat %>% mutate(pair_type = factor(pair_type, levels = c('Intra-facility pair', 'Inter-facility pair'))) %>%
  ggplot(aes(x = reorder(loc1, -n, median), y = n, fill = pair_type)) +
  geom_col(position = 'dodge') +
       scale_fill_manual(values=c("cadetblue", "salmon")) + 
         labs(y = "Number of closely related \nisolate pairs", x = "Facility", fill = '', color = '') +
  facet_grid(factor(snvthresh, levels = c('< 7 SNVs','< 11 SNVs'))~.)

## ---- message=FALSE-----------------------------------------------------------
# get dates 
dates <- metadata %>% select(sample_id, collection_date) %>% deframe()

time1 <- snv_dists %>% 
  mutate(date1 = dates[sample1], date2 = dates[sample2]) %>% 
  filter(date1 == date2) %>% 
  ggplot(aes(x = pairwise_dist, fill = pair_type)) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 30) +
  labs(x = "Pairwise SNV distance", y = "Count", fill = "") +
  facet_grid(~date1) + 
  theme(legend.justification = "top", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

time2 <- snv_dists %>% 
  mutate(date1 = dates[sample1], date2 = dates[sample2]) %>% 
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

## ---- message=FALSE-----------------------------------------------------------
# get summary of pairwise SNV distance data
snv_summary <- summarize_pairs(snv_dists,summary_fns = c('min'), threshs = c(7, 11)) %>% filter(loc1 != loc2)
snv_summary

## -----------------------------------------------------------------------------
head(pt_trans_df)

## -----------------------------------------------------------------------------
pt_trans <- get_patient_flow(edge_df = pt_trans_df)
head(pt_trans)

## -----------------------------------------------------------------------------
# have to convert fsp matrix to long form first
fsp_long <- make_long_form(fsp)

pair_info <- merge_inter_summaries(pt_trans, snv_summary, fsp_long) 
pair_info

#%>% ggplot(mapping = aes(x = sum_transfers, y = under_7)) + geom_jitter(position = position_jitter(width = 0.2, height = 0.1), alpha = 0.2, color = "blue") + labs(y = "Number of Closely Related sample Pairs", x = "Number of Transfers", title = "Number of Patient Transfers between Facilities \nvs. Number of Closely Related sample Pairs at Facilities \nfor Each Facility Pair") + geom_smooth(method = "lm")



## -----------------------------------------------------------------------------
pair_info %>% 
  ggplot(aes(x=sum_pt_trans_metric,y=fsp)) +
  geom_point() + geom_smooth(method='lm') + scale_x_log10() +
  labs(x = 'Patient flow', y = 'Fsp') 

## -----------------------------------------------------------------------------
pair_info %>% 
  ggplot(aes(x=sum_pt_trans_metric,y=under_11)) +
  geom_point() + geom_smooth(method='lm') + scale_x_log10() +
  labs(x = 'Patient flow', y = '# closely related\npairs (< 11 SNVs)') 

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
snv_geo_thresh <- snv_geo %>% filter(under_7 > 0) 

facil_geo %>% 
  ggplot(aes(x=long, y=lat)) +
  geom_curve(aes(x = facil_1_long, y = facil_1_lat, xend = facil_2_long, yend = facil_2_lat,
                 alpha = under_7),
             data = snv_geo_thresh, curvature = 0.33) +
  geom_point(aes(size=n_isolates), alpha=1, color = "lightblue") + 
  theme_map() + 
  labs(size = 'Number of samples', alpha = 'Number of closely related\npairs (< 7 SNVs)')

## -----------------------------------------------------------------------------
snv_geo <- snv_geo %>% 
  mutate(geo_dist = sqrt((facil_1_long-facil_2_long)^2+(facil_1_lat-facil_2_lat)^2))

## -----------------------------------------------------------------------------
snv_geo %>% 
  ggplot(aes(x = geo_dist, y = under_11)) + geom_point() + geom_smooth(method = 'lm') +
  labs(x = 'Physical distance between facilities', y = '# closely related pairs\n(< 11 SNVs)')

## -----------------------------------------------------------------------------
# merge snv_dists and metadata, and find time difference
snv_and_metadat <- left_join(snv_dists, metadata, by = c('sample1' = 'sample_id')) %>% 
  rename(patient_id_1 = patient_id, collection_date_1 = collection_date) %>% 
  left_join(metadata, by = c('sample2' = 'sample_id')) 
head(snv_and_metadat)

## ---- eval = FALSE------------------------------------------------------------
#  doFuture::registerDoFuture()
#  future::plan(future::multisession, workers = 4)
#  
#  clusters <- get_clusters(tr, locs)

