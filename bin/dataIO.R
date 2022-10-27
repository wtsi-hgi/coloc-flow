#!/usr/bin/env Rscript
library(data.table)
requireNamespace('dplyr')
requireNamespace('tidyr')

load_GWAS <- function(GWAS){
  message(paste('Reading GWAS:', GWAS))

  GWAS_ext = tools::file_ext(GWAS)
  GWAS_name = tools::file_path_sans_ext(basename(GWAS))

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

    'REF' = "effect_allele",
    'Allele1' = "effect_allele",

    'ALT' = "other_allele",
    'Allele2' = "other_allele",

    'CHROM' = "chromosome",
    'CHR' = "chromosome",

    'ID' = "variant_id",
    'MarkerName' = "variant_id",

    'POS' = "base_pair_location",
    'pos' = "base_pair_location",

    'BETA' = "beta",
    'Effect' = "beta",

    'SE' = "standard_error",
    'StdErr' = "standard_error",

    # 'NStudies' = "N",
    'OBS_CT' = "N"
  )

  #Gwas col rename
  table_rules <- renaming_rules[intersect(names(renaming_rules), colnames(map))]
  map <- dplyr::rename(map, !!!setNames(names(table_rules), nm=table_rules))

  return_list <- list("map" = map, "GWAS_name" = GWAS_name)
  return(return_list)

}

#### Eqtl data
load_eqtl <- function(eqtl.file, marker.file, build = 'hg38'){
    message(paste("Reading eQTL:", eqtl.file))
    build_choices <- c('hg19', 'hg38')
    build <- match.arg(build, choices = build_choices)

    # https://zenodo.org/record/6104982
    eqtl <- fread(eqtl.file, col.names = c("gene", "SNP", "distance_to_TSS", "p", "beta"))
    eqtl$se <- abs(eqtl$beta) / sqrt(qchisq(eqtl$p, df = 1, lower.tail = F))

    position_column_name <- paste0('SNP_id_', build)
    another_column_name <- paste0('SNP_id_', setdiff(build_choices, build))

    snps <- fread(marker.file, drop = another_column_name)

    single_eqtl <- dplyr::inner_join(eqtl, snps, by = 'SNP')
    single_eqtl <- tidyr::separate(single_eqtl, col = position_column_name, sep = ':', convert = T,
                                   into = c('chromosome', 'base_pair_location'))
    single_eqtl$chromosome <- as.integer(gsub('chr', '', single_eqtl$chromosome))

    #eQTL col rename
    single_eqtl <- dplyr::rename(single_eqtl,
      p_value = p,
      variant_id = SNP,
      standard_error = se
    )
    return (single_eqtl)
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
