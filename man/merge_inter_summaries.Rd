% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_inter_pairs.R
\name{merge_inter_summaries}
\alias{merge_inter_summaries}
\title{Merge summarized data about facility pairs}
\usage{
merge_inter_summaries(
  patient_flow = NULL,
  inter_pair_summary = NULL,
  fsp_long = NULL
)
}
\arguments{
\item{patient_flow}{output of get_patient_flow function}

\item{inter_pair_summary}{output of summarize_pairs function}

\item{fsp}{output of get_facility_fsp function}
}
\value{
merged dataframe
}
\description{
Merge summarized data about facility pairs
}
\examples{
locs <- metadata \%>\% dplyr::select(isolate_id, facility) \%>\% tibble::deframe()
snv_dists <- get_snv_dists(dists, locs)
inter_pair_summary <- summarize_pairs(snv_dists)
patient_flow <- get_patient_flow(edge_df = pt_trans_df)
fsp_long <- make_long_form(fsp)
merge_inter_summaries(patient_flow, inter_pair_summary, fsp_long)
}
