# regentrans: R package for investigating regional pathogen transmission using genomic data

<img src="man/figures/regentrans_logo.png"  width = "200" />

Regentrans is a work in progress. 

It can be installed using the command `devtools::install_github('Snitkin-Lab-Umich/regentrans')` in R 

An introductory vignette can be found [here](articles/Introduction.html).

Questions regentrans can help investigate:

| Question | Method | regentrans function(s) | Required Data | Optional Data
|---|---|---|---|---|
| Is transmission occurring in the region of interest? | Pairwise SNV distances | `get_snv_dists()` | A pairwise SNV distance matrix (can be created using `ape::dist.dna()` on fasta file of variants), isolate location information | Isolate patient information, patient transfer network |
| Is transmission occurring within and/or between locations? | Phylogenetic clustering of isolates from the same location | `get_clusters()` | Phylogenetic tree, isolate location information | |
| | Pairwise SNV distances within and between facilities | `get_snv_dists()` | A pairwise SNV distance matrix (can be created using `ape::dist.dna()` on fasta file of variants), isolate location information | Isolate patient information, patient transfer network |
| | Fraction of intra-facility pairs for different pairwise SNV distance thresholds | `get_frac_intra()` | A pairwise SNV distance matrix (can be created using `ape::dist.dna()` on fasta file of variants), isolate location information | Isolate patient information, patient transfer network |
| What locations is transmission occurring within/between? | Closely related pairs within and between facilities | `get_snv_dists()` | A pairwise SNV distance matrix (can be created using `ape::dist.dna()` on fasta file of variants), isolate location information | Isolate patient information, patient transfer network |
| | Population-level similarity between locations | `get_facility_fsp()` | Fasta file of variants, isolate location information | |
| Have transmission dynamics changed over time? | Methods above but split over time | | | |
| Is transmission occurring along paths of higher patient/person flow? | Compare patient/person flow between locations to inter-location pairwise SNV distances | `patient_transfer()`| Patient transfer network, a pairwise SNV distance matrix (can be created using `ape::dist.dna()` on fasta file of variants), isolate location information | Isolate patient information |
| | Compare patient/person flow between locations to populatopn-level similarity between locations | `patient_transfer()`, `get_facility_fsp()`| Patient transfer network, fasta file of variants, isolate location information | |
| Are there any observable geographic trends in prevalence/transmission? | Visualize prevalence and closely related pairs over time | Plot code in vignette, `patient_transfer()` *or* `get_facility_fsp()` | Geographic locations of each facility, patient transfer network, fasta file of variants, isolate location information | |
