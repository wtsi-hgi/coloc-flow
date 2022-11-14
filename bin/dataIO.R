#!/usr/bin/env Rscript
library(data.table)
requireNamespace('dplyr')
requireNamespace('tidyr')

load_GWAS <- function(GWAS){
  message(paste('Reading GWAS:', GWAS))

  GWAS_ext = tools::file_ext(GWAS)
  GWAS_name = tools::file_path_sans_ext(basename(GWAS), compression = T)

  if(GWAS_ext=='zip'){
    print('zip')
    map = fread(paste0('unzip -p ', GWAS))

  } else if (GWAS_ext=='gz'){
    map = fread(paste0('gunzip -cq ', GWAS))
  } else{
    map = fread(GWAS)
  }

  renaming_rules <- list(
    'P' = "p_value",
    'P-value' = "p_value",

    'REF' = "reference_allele",
    'ALT' = "alternative_allele",

    'Allele1' = "effect_allele",
    'A1' = "effect_allele",
    'A_EFF' = "effect_allele",

    'Allele2' = "other_allele",
    'AX' = 'other_allele',
    'A_NONEFF' = 'other_allele',

    'CHROM' = "chromosome",
    'CHR' = "chromosome",

    'ID' = "variant_id",
    'MarkerName' = "variant_id",
    'SNP' = "variant_id",

    'POS' = "base_pair_location",
    'pos' = "base_pair_location",
    'BP' = "base_pair_location",

    'BETA' = "beta",
    'Effect' = "beta",

    'SE' = "standard_error",
    'StdErr' = "standard_error",

    'Freq1' = 'eaf',
    'A1_FREQ' = 'eaf',
    'Freq_EFF' = 'eaf',

    'OBS_CT' = "N",
    'N' = 'N'
  )

  col.names <- unique(renaming_rules)

  # Gwas col rename
  table_cols <- c(
    setNames(names(renaming_rules), nm=renaming_rules),
    unlist(col.names)
  )
  map <- dplyr::select(map, !!!dplyr::any_of(table_cols))

  return_list <- list("map" = map, "GWAS_name" = GWAS_name)
  return(return_list)
}

#### Eqtl data
load_eqtl <- function(eqtl.file, marker.file, marker.data, build = 'hg38'){
    message(paste("Reading eQTL:", eqtl.file))

    # https://zenodo.org/record/6104982
    eqtl <- fread(eqtl.file, col.names = c("gene", "SNP", "distance_to_TSS", "p", "beta"))
    eqtl$se <- abs(eqtl$beta) / sqrt(qchisq(eqtl$p, df = 1, lower.tail = F))

    if(missing(marker.data)){
      marker.data <- read_eqtl_marker_file(marker.file, build)
    }

    single_eqtl <- dplyr::inner_join(eqtl, marker.data, by = 'SNP')

    #eQTL col rename
    single_eqtl <- dplyr::rename(single_eqtl,
      p_value = p,
      variant_id = SNP,
      standard_error = se
    )
    return (single_eqtl)
}

read_eqtl_marker_file <- function (marker.file, build = 'hg38'){
    message(paste("Reading file:", marker.file))
    build_choices <- c('hg19', 'hg38')
    build <- match.arg(build, choices = build_choices)

    position_column_name <- paste0('SNP_id_', build)
    another_column_name <- paste0('SNP_id_', setdiff(build_choices, build))

    snps <- fread(marker.file, drop = another_column_name)
    single_eqtl <- tidyr::separate(snps, col = position_column_name, sep = ':', convert = T,
                                   into = c('chromosome', 'base_pair_location'))
    single_eqtl$chromosome <- as.integer(gsub('chr', '', single_eqtl$chromosome))

    return(single_eqtl)
}

# reads plink's frq(x) file
read_freqs <- function (filename){
    message(paste("Reading freqs: ", filename))
    ext <- tools::file_ext(filename)
    stopifnot(ext %in% c('frq', 'frqx'))

    freqs <- fread(filename)
    if(ext == 'frq'){
        freqs <- dplyr::rename(freqs, FreqA1 = MAF, N = NCHROBS)
        freqs$N <- freqs$N/2L
    }

    if(ext == 'frqx'){
        freqs$N <- rowSums(freqs[, c("C(HOM A1)", "C(HET)", "C(HOM A2)")])
        freqs$FreqA1 <- (freqs$'C(HOM A1)'*2 + freqs$'C(HET)') / (2*freqs$N)
    }

    return(freqs)
}

# if gwas data is missing eaf, this function will add it from freq table
add_freq <- function (df, freq){
    require(dplyr)

    m <- inner_join(
      df %>% mutate(chromosome = as.integer(as.character(chromosome))),
      freq %>% select(-OBS_CT),
      by = c('chromosome'='#CHROM', 'variant_id'='ID')
    ) %>%
      rowwise() %>%
      filter(all(c(effect_allele, other_allele) %in% c(REF, ALT))) %>%
      ungroup()

    df_freq <- mutate(m, eaf = ifelse(ALT == effect_allele, ALT_FREQS, 1-ALT_FREQS)) %>%
      select(-REF, -ALT, -ALT_FREQS)

    return(df_freq)
}

add_position <- function (df, variant_colname = 'variant_id', genotype_prefix){
    stopifnot(variant_colname %in% colnames(df))
    bim_filename <- paste(genotype_prefix, 'bim', sep='.')
    bim <- fread(bim_filename, col.names = c('chromosome', 'rsid', 'base_pair_location'), select = c(1, 2, 4))
    m <- dplyr::inner_join(df, bim, by = setNames('rsid', variant_colname))
    return(m)
}
