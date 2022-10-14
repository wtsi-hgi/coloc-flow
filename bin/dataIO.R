library(data.table)
requireNamespace('dplyr')

load_GWAS <- function(GWAS){

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

load_eqtl <- function(eqtl.file, Full_GWAS_Sum_Stats){
  #### Eqtl data
  eqtl = fread(eqtl.file, col.names = c("gene","SNP","TSS_dist","p","beta"))
  eqtl$se = abs(eqtl$beta)/sqrt(qchisq(eqtl$p,df = 1,lower.tail = F))
  Celltype = basename(eqtl.file)
  Celltype = gsub('\\.','_',Celltype)
  single_eqtl1 = eqtl

  # Here we diver a bit and add the positional info to the SNPs if available.
  snp_matches <- match(single_eqtl1$SNP,Full_GWAS_Sum_Stats$variant_id)
  single_eqtl1$ea_allele = Full_GWAS_Sum_Stats$effect_allele[snp_matches]
  single_eqtl1$oth_allele = Full_GWAS_Sum_Stats$other_allele[snp_matches]
  single_eqtl1$chromosome = Full_GWAS_Sum_Stats$chromosome[snp_matches]
  single_eqtl1$base_pair_location = Full_GWAS_Sum_Stats$base_pair_location[snp_matches]
  single_eqtl1 = single_eqtl1[!is.na(single_eqtl1$ea_allele)]

  #eQTL col rename
  single_eqtl1 <- dplyr::rename(single_eqtl1,
      p_value = p,
      effect_allele = ea_allele,
      other_allele = oth_allele,
      variant_id = SNP,
      standard_error = se
  )
  return (single_eqtl1)
}

# reads plink's frq(x) file
read_freqs <- function (filename){
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
