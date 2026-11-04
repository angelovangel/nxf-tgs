#!/usr/bin/env Rscript

# combine wf-bacterial-genomes results.json and samplesheet_validated.csv
# produce sample_status.csv summary for every user (same base columns as sample_status.R)

# arg[1] is user
# arg[2] is results.json (from wf-bacterial-genomes out_dir)
# arg[3] is samplesheet_validated
#
# samples missing from results.json, with sample_pass false or without contigs are 'fail'
library(dplyr)
library(vroom)
library(jsonlite)

arg <- commandArgs(trailingOnly = T)

res <- fromJSON(arg[2], simplifyVector = F)

# NULL/empty safe accessor
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# top sourmash species call with its ANI, e.g. 'Acinetobacter baumannii (98.7%)'
species_with_ani <- function(matches) {
  if (length(matches) == 0) return(NA_character_)
  species <- matches[[1]]$species_call %||% NA_character_
  ani <- matches[[1]]$query_containment_ani
  if (is.na(species) || is.null(ani)) return(species)
  ani <- as.numeric(ani)
  if (ani <= 1) ani <- ani * 100 # sourmash reports a fraction
  sprintf('%s (%.1f%%)', species, ani)
}

# other species (different genus than the top call) explaining >= 5% of the assembly - mixed sample
# e.g. 'Acinetobacter baumannii (21%)'; related species of the same genus at ~1-2% are normal
other_species <- function(matches, min_fraction = 0.05) {
  if (length(matches) < 2) return(NA_character_)
  genus <- function(m) sub(' .*', '', m$species_call %||% '')
  top_genus <- genus(matches[[1]])
  other <- Filter(function(m) genus(m) != top_genus && (m$f_unique_to_query %||% 0) >= min_fraction, matches[-1])
  if (length(other) == 0) return(NA_character_)
  paste(vapply(other, function(m) sprintf('%s (%.0f%%)', m$species_call, 100 * m$f_unique_to_query), character(1)), collapse = '; ')
}

# ResFinder classes of the drugs predicted resistant, carbapenems reported separately from other beta-lactams
resistance_classes <- function(amr) {
  resistant <- Filter(function(p) isTRUE(p$amr_resistant), amr$phenotypes %||% list())
  if (length(resistant) == 0) return(NA_character_)
  classes <- unique(unlist(lapply(resistant, function(p) p$amr_classes)))
  classes <- setdiff(classes, 'under_development')
  drugs <- vapply(resistant, function(p) p$drug %||% '', character(1))
  if (any(drugs %in% c('imipenem', 'meropenem', 'ertapenem'))) classes <- c('carbapenem', setdiff(classes, 'beta-lactam'))
  paste(sort(unique(classes)), collapse = ', ')
}

# json with one entry per gene for the expandable details in the user report: associated drugs of all hits,
# identity and coverage of the best hit (highest ResFinder grade) and the contigs with a hit
amr_gene_details <- function(regions) {
  if (length(regions) == 0) return(NA_character_)
  hits <- bind_rows(lapply(regions, function(g) data.frame(
    gene = g$gene_name %||% '',
    drug = unlist(g$phenotypes) %||% NA_character_,
    grade = g$grade %||% NA_real_,
    identity = g$identity %||% NA_real_,
    coverage = g$coverage %||% NA_real_,
    contig = g$query_contig %||% NA_character_,
    stringsAsFactors = F
  )))
  details <- hits %>%
    arrange(desc(grade), desc(identity), desc(coverage)) %>%
    group_by(gene) %>%
    summarise(
      drugs = paste(sort(unique(na.omit(drug))), collapse = ', '),
      identity = round(first(identity), 1),
      coverage = round(first(coverage), 1),
      contigs = paste(unique(na.omit(contig)), collapse = ', ')
    )
  as.character(toJSON(details, dataframe = 'rows'))
}

# 'ST19 (abaumannii_2)', '-' when no ST could be assigned
mlst <- function(st) {
  type <- st$sequence_type %||% NA_character_
  if (is.na(type) || type == '-') return(type)
  sprintf('ST%s (%s)', type, st$detected_species %||% '?')
}

parse_sample <- function(s) {
  r <- s$results
  contigs <- r$assembly$contig %||% list()
  flye <- r$assembly$flye_stats %||% list()
  lengths <- vapply(contigs, function(x) as.numeric(x$length %||% NA_real_), numeric(1))
  depths <- vapply(contigs, function(x) as.numeric(x$coverage$mean %||% NA_real_), numeric(1))
  n_contigs <- length(contigs)
  matches <- r$species_identification$detected_matches %||% list()
  regions <- r$antimicrobial_resistance$seq_regions %||% list()
  data.frame(
    sample = as.character(s$alias),
    n_seqs = r$fastq$n_seqs %||% NA_real_,
    pass_fail = if (isTRUE(s$sample_pass) && n_contigs > 0) 'pass' else 'fail',
    assembly_size = if (n_contigs > 0) sum(lengths) else NA_real_,
    n_contigs = n_contigs,
    n_circular = sum(vapply(flye, function(x) identical(x$circular, 'Y'), logical(1))),
    largest_contig = if (n_contigs > 0) max(lengths) else NA_real_,
    # length-weighted mean depth over all contigs
    mean_depth = if (n_contigs > 0) round(sum(lengths * depths) / sum(lengths)) else NA_real_,
    n_plasmids = r$plasmid_identification$n_plasmid_contigs %||% NA_integer_,
    species = species_with_ani(matches),
    other_species = other_species(matches),
    mlst_st = mlst(r$sequence_typing),
    resistance = resistance_classes(r$antimicrobial_resistance),
    # distinct genes, the same gene can be hit on several contigs
    n_amr_genes = length(unique(vapply(regions, function(g) g$gene_name %||% '', character(1)))),
    amr_genes = amr_gene_details(regions),
    stringsAsFactors = F
  )
}

extra_cols <- c('n_contigs', 'n_circular', 'largest_contig', 'mean_depth', 'n_plasmids', 'species', 'other_species', 'mlst_st', 'resistance', 'n_amr_genes', 'amr_genes')

df1 <- bind_rows(lapply(res$samples, parse_sample))
if (nrow(df1) == 0) {
  df1 <- data.frame(sample = character(), n_seqs = numeric(), pass_fail = character(), assembly_size = numeric())
  for (col in extra_cols) df1[[col]] <- logical()
}

#samplesheet_validated
df2 <- vroom(
  arg[3],
  col_names = T,
  trim_ws = T, delim = ","
  ) %>%
  dplyr::filter(user == arg[1]) %>%
  dplyr::mutate(sample = as.character(sample))

# get_maxbin.sh (obs_size, nreads) is not run for wf-bacterial-genomes, add the missing columns
for (col in c('dna_size', 'obs_size', 'nreads')) {
  if (!col %in% names(df2)) df2[[col]] <- NA_real_
}

df2 %>%
  left_join(df1, by = join_by(sample)) %>%
  dplyr::mutate(
    pass_fail = coalesce(pass_fail, 'fail'),
    nreads = coalesce(nreads, n_seqs), # reads after wf-bacterial-genomes length filtering
    assembly_quality = NA_real_ # no per-sample assembly quality in wf-bacterial-genomes
  ) %>%
  dplyr::select(
    c('user', 'sample', 'barcode', 'validate', status = 'pass_fail', 'nreads', user_size = 'dna_size', 'obs_size', 'assembly_size', 'assembly_quality'),
    all_of(extra_cols)
  ) %>%
  write.csv(file = 'sample-status.csv', row.names = F)
