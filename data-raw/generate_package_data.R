#subset Penn dataset for github
library(ape)

# metadata
kp_metadata <- read.delim(system.file("data-raw", "kp_metadata.csv", package = "regentrans"), sep = ",")

# days between isolates
days_between_isolates <- read.delim(system.file("data-raw", "days_between_isolates.csv", package = "regentrans"), sep = ",")

# alignment
kp_aln <- ape::read.dna(system.file("data-raw", "kp.fasta", package = "regentrans"),
                       format = "fasta")

# tree
kp_tr <- ape::read.tree(system.file("data-raw", "kp.treefile", package = "regentrans"))

# get pairwise SNV distances
kp_dists <- ape::dist.dna(x = kp_aln, as.matrix = TRUE, model = "N", pairwise.deletion = TRUE) %>% data.frame()

# get STs
st258_ids <- read.delim(system.file("data-raw", "kp_st.csv", package = "regentrans"), sep = ",") %>%
  dplyr::filter(st == 'ST258') %>% dplyr::select(isolate_id) %>% tibble::deframe()
# include only ones in alignment
st258_ids <- st258_ids[st258_ids %in% rownames(kp_aln)]

# subset to ST258
metadata <- kp_metadata %>% dplyr::filter(isolate_id %in% st258_ids)
aln <- kp_aln[st258_ids,]
tr <- ape::keep.tip(kp_tr, st258_ids)
dists <- kp_dists[st258_ids, st258_ids]
days_between_isolates <- days_between_isolates %>%
  dplyr::filter(isolate1 %in% st258_ids & isolate2 %in% st258_ids)

# create mini dataset
metadata_mini <- metadata[1:10,]
aln_mini <- aln[metadata_mini$sample_id,]
tr_mini <- ape::keep.tip(tr, metadata_mini$sample_id)
dists_mini <- dists[metadata_mini$sample_id, metadata_mini$sample_id]

# save data to package
usethis::use_data(metadata, overwrite = TRUE)
usethis::use_data(days_between_isolates, overwrite = TRUE)
usethis::use_data(aln, overwrite = TRUE)
usethis::use_data(tr, overwrite = TRUE)
usethis::use_data(dists, overwrite = TRUE)

usethis::use_data(metadata_mini, overwrite = TRUE)
usethis::use_data(aln_mini, overwrite = TRUE)
usethis::use_data(tr_mini, overwrite = TRUE)
usethis::use_data(dists_mini, overwrite = TRUE)

# include example fsp (because it takes a while to run)
# run these 2 lines if need to update fsp
# fsp <- get_facility_fsp(aln, locs)
# write.csv(fsp, 'data-raw/kp_fsp.csv', quote = FALSE)
fsp <- read.delim(system.file("data-raw", "kp_fsp.csv", package = "regentrans"), sep = ",", row.names = 1)

# save fsp to package
usethis::use_data(fsp, overwrite = TRUE)

# include true indirect patient flow network (because can't share raw data)
pt_flow <- structure(list(D = c(NA, 1.71122994652408e-05, 1.40967224149537e-05,
                                0.000474074074074073, 0.000181771667856037, 5.96955526813254e-05,
                                6.25341182920926e-05, 0.000329896907216494, 1.79559003088417e-05,
                                0.00424361990053091, 3.76915038187176e-05),
                          H = c(0.000260310849248276,NA, 2.11840865735182e-05, 0.00049223929181057,
                                0.00264038021475091, 2.78373043591672e-05, 0.00195723974827721,
                                6.31601199268895e-05, 2.71526055905228e-05, 0.000254772320540862, 0.000117836481565193
                                ),
                          G = c(5.06444506343215e-06, 1.02719860784377e-05, NA, 2.16590859865714e-05,
                                4.2735042735043e-05, 0.000935486692948707, 0.00139281644367536,
                                7.57161860420929e-05, 0.00537240537240539, 1.84891380716003e-05,
                                0.00203488372093025),
                          C = c(0.00266006237387636, 0.000655364660370925, 0.000188781104829633,
                                NA, 6.54568436258582e-05, 0.000831758034026462, 0.00270305876961852,
                                0.00163163643144462, 2.37908450468836e-05, 0.0106382978723404, 0.00105009489561482),
                          T = c(8.28622018387643e-06, .00279720279720279, 2.82027726589947e-06, 8.98957209636826e-05,
                                NA, 2.04140231967224e-05, 0.000378469301934397, 3.00762475842331e-06,
                                2.1399518254864e-06, 5.4190701474378e-06, 5.16956162117448e-05),
                          F = c(4.58631443771784e-05, 3.26601235133297e-05, 0.00218315188915409,
                                0.000562962962962956, 4.35085276714237e-05, NA, 9.54505264850083e-05,
                                0.00705090540035253, 0.00098948670377241, 4.48873327946856e-05,
                                0.000102291325695581),
                          I = c(1.24081794719078e-05, 0.000768716577540108, 0.000142391135504582, 0.00320593942461823,
                                0.00243902439024388, 0.000596273291925468, NA, 1.96367206676486e-05,
                                1.82868846463318e-05, 3.41057385597683e-05, 0.0005686517783292),
                          E = c(0.000504494588148968, 0.0001241002730206, 0.00143820224719103, 0.00171765968867417,
                                4.43880103948647e-05, 0.0111304347826088, 0.000190901052970016,
                                NA, 0.000115440115440116, 0.000493760660741536, 0.000234603040455405),
                          P = c(7.93788907105321e-06, 7.61436816906775e-06, 0.0104238932104101, 7.76601454985673e-06,
                                2.59839246119733e-05, 0.000405496733498536, 0.000333235091855784, 3.00854426571462e-05,
                                NA, 0.000172047944027069, 0.00287122366962305),
                          M = c(0.00329188929123093, 0.000318606627017839, 4.48996627051994e-05, 0.00104341827275471,
                                4.35085276714237e-05,  9.83671060397392e-06, 9.75191137462946e-06,
                                7.67266994963933e-05, 3.97355203763749e-05, NA, 0.000249754373797667),
                          O = c(6.22642852060063e-06, 7.51936235807199e-06, 0.000125335456663423, 6.12688784731793e-06,
                                3.18126868995355e-05, 6.00541032875977e-05, 1.09443806574985e-05,
                                1.83618814685465e-05, 0.00108745087419629, 3.94496051301495e-06, NA)),
                     class = "data.frame", row.names = c("D", "H", "G", "C", "T", "F", "I", "E", "P", "M", "O"))
diag(pt_flow) <- NA

# save patient flow to package
usethis::use_data(pt_flow, overwrite = TRUE)

#made up patient flow as example (using real data as template)
pt_trans_df <- round(pt_flow * 1e6) %>%
  make_long_form(col_names = c('source_facil','dest_facil','n_transfers')) %>%
  dplyr::filter(!is.na(n_transfers))

# save made up patient transfer dataframe to package
usethis::use_data(pt_trans_df, overwrite = TRUE)

# deidentified longitude and latitude
# these are deidentified
longs <- c(A = NA, B = NA, C = -930605.082601, D = -930605.009738, E = -930604.936812,
           F = -930604.912683, G = -930604.930899, H = -930605.236274, I = -930605.065595,
           J = NA, K = NA, L = NA, M = -930605.109771, N = NA, O = -930605.187588,
           P = -930604.896185, Q = NA, R = NA, S = NA, T = -930605.237532,
           U = NA)
lats <- c(A = NA, B = NA, C = -740730.013261, D = -740730.378523, E = -740729.967478,
          F = -740729.917235, G = -740729.643727, H = -740729.995281, I = -740729.894068,
          J = NA, K = NA, L = NA, M = -740730.297362, N = NA, O = -740729.22271,
          P = -740729.566327, Q = NA, R = NA, S = NA, T = -740729.891502,
          U = NA)

# make longitude/latitude data frame
facil_coord <- bind_cols(facil=names(longs), long=longs, lat=lats) %>%
  filter(!is.na(long) & !is.na(lat))

# save deidentified longitude and latitude to package
usethis::use_data(facil_coord, overwrite = TRUE)
