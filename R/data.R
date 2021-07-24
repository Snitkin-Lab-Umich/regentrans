#' Klebiella pneumoniae ST258 isolate metadata
#'
#' Metadata including isoalte ID, patient ID, facility, and isolate collection
#' date of Klebsiella pneumoniae ST258 isolates collected from long-term acute
#' care hospitals across the US.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A data frame with 413 rows and 4 variables.
"metadata"

#' Klebiella pneumoniae ST258 isolate variant alignment
#'
#' Variant alignment for Klebsiella pneumoniae ST258 isolates in the study.
#' samples were collected from long-term acute care hospitals across the US.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A DNAbin object with 413 rows and 83976 nucelotide sites.
"aln"

#' Klebiella pneumoniae ST258 isolate phylogenetic tree
#'
#' Phylogenetic tree of Klebsiella pneumoniae ST258 isolates in the study.
#' samples were collected from long-term acute care hospitals across the US.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A phylo object with 413 tips.
"tr"

#' Klebiella pneumoniae ST258 isolate pairwise SNV distance
#'
#' Pairwise single nucleotide variant (SNV) distance matrix for
#' Klebsiella pneumoniae ST258 isolates in the study.
#' samples were collected from long-term acute care hospitals across the US.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A data frame with 413 rows and 413 columns.
"dists"

#' Klebiella pneumoniae ST258 gene flow
#'
#' Gene flow (fsp) between facilities for Klebsiella pneumoniae ST258 isolates.
#' samples were collected from long-term acute care hospitals across the US.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A data frame with 413 rows and 413 columns.
"fsp"

#' Patient flow between facilities
#'
#' Patient flow between California facilities calculated from a patient
#' transfer network of the number of patients transferred between
#' facilities over the course of a year.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A data frame with 11 rows and 11 columns.
"pt_flow"

#' Made up patient transfer dataframe between facilities
#'
#' Made up patient transfers between facilities based on the true patient
#' flow network. We can't provide the original data here.
#' Corresponding manuscript: https://aac.asm.org/content/63/11/e01622-19
#'
#' @format A data frame with 11 rows and 11 columns.
"pt_trans_df"

#' Subsetted metadata
#'
#' Small version for examples.
#'
#' @format A data frame
"metadata_mini"

#' Subsetted metadata
#'
#' Subset of `metadata`. Small version for examples.
#'
#' @format A data frame.
"metadata_mini"

#' Subsetted alignment
#'
#' Subset of `aln`. Small version for examples.
#'
#' @format A DNAbin object.
"aln_mini"

#' Subsetted tree
#'
#' Subset of `tr`. Small version for examples.
#'
#' @format A phylo object.
"tr_mini"

#' Subsetted pairwise SNV distances
#'
#' Subset of `dists`. Small version for examples.
#'
#' @format A data frame.
"dists_mini"


#' Deidentified facility longitude and latitude
#'
#' Deidentified facility longitude and latitude for plotting
#'
#' @format A data frame.
"facil_coord"

