#!/usr/bin/env Rscript
library(data.table)
requireNamespace('dplyr')
requireNamespace('tidyr')
library("stringr") 
library(tidyr)


#### Eqtl data
load_eqtl <- function(eqtl.file, marker.file, marker.data, build = 'hg38',eqtl_significance_threshold=10){
    message(paste("Reading file:", eqtl.file))
    
    # https://zenodo.org/record/6104982
    # eqtl <- fread('/scratch/cellfunc/shared/HUVEC_RNAseq/eQTLs_norm_counts/TensorQTL_eQTLS/general/nom_output/cis_nominal1.cis_qtl_pairs.chr1.tsv')
    # eqtl.file=eQTL
    # marker.file = eqtl_marker_file
    eqtl <- fread(eqtl.file)
    #eqtl= map
    # check if file has headers, if it doesnt then ve assume the order is as per:
    # "gene", "SNP", "distance_to_TSS", "p", "beta"
    if ('V1' %in% names(eqtl)){
        names(eqtl)=c("gene", "SNP", "distance_to_TSS", "p", "beta")
        eqtl$se <- abs(eqtl$beta) / sqrt(qchisq(eqtl$p, df = 1, lower.tail = F))
    }else{
      renaming_rules <- list(
        'P' = "p",
        'P-value' = "p",
        'pval_nominal' = "p",
        'tss_distance' = "distance_to_TSS",
        'variant_id' = "SNP",
        'phenotype_id' = "gene",
        'slope' = 'beta',
        'slope_se'='se',


        'REF' = "reference_allele",
        'ALT' = "alternative_allele",

        'Allele1' = "effect_allele",
        'A1' = "effect_allele",
        'A_EFF' = "effect_allele",
        'EA' = "effect_allele",

        'Allele2' = "other_allele",
        'AX' = 'other_allele',
        'A_NONEFF' = 'other_allele',
        'NEA' = 'other_allele',

        'CHROM' = "chromosome",
        'CHR' = "chromosome",

        'ID' = "SNP",
        'MarkerName' = "SNP",
        'rsID' = "SNP",

        'POS' = "base_pair_location",
        'pos' = "base_pair_location",
        'BP' = "base_pair_location",

        'BETA' = "beta",
        'Effect' = "beta",

        'SE' = "se",
        'StdErr' = "se",

        'Freq1' = 'eaf',
        'A1_FREQ' = 'eaf',
        'Freq_EFF' = 'eaf',
        'EAF' = 'eaf',

        'OBS_CT' = "N",
        'N' = 'N'

      )

      col.names <- unique(renaming_rules)

      # Gwas col rename
      table_cols <- c(
        setNames(names(renaming_rules), nm=renaming_rules),
        unlist(col.names)
      )
      eqtl <- dplyr::select(eqtl, !!!dplyr::any_of(table_cols))

      # Check if beta exists, if not calculate this.
      if (!'beta' %in% names(eqtl)){
        eqtl$beta = eqtl$se * sqrt(qchisq(eqtl$p, df = 1, lower.tail = F))
      }
    }

    if(missing(marker.data)){
      marker.data <- read_eqtl_marker_file(marker.file, build)
    }

    # Check if here we have a variant id in position format: chr1:68866536:T:G
    #  if so convert this to rsid.
    # eqtl$variant_id[0]
    
    # For the eqtls have to add chromosome and base_pair_location for the variants.
    # We only want to add the positions in cases when there is no existing positions already available.
    if (!("chromosome" %in% colnames(eqtl)) & !('base_pair_location' %in% colnames(eqtl))){
        single_eqtl <- dplyr::inner_join(eqtl, marker.data, by = 'SNP')
    } else{
       single_eqtl <- eqtl
    }
    

    if(build != 'hg38'){
        message(paste('Convert positions from', build, 'to hg38'))
        ch <- load_chain_file(from = build, to = 'hg38')
        single_eqtl <- lift_over_df(single_eqtl, chain = ch)

    }

    # single_eqtl = single_eqtl[single_eqtl$p < eqtl_significance_threshold]

    if (sum(str_detect(single_eqtl$SNP, ':')) > 0){
        # This part checks for the ids that needs to be converted and convers them to rsids where available.
        # quite often there are no rsids associated. 
        # For this we could consider converting GWAS loci to chr positons.
        to_fix = single_eqtl[str_detect(single_eqtl$SNP, ':')]
        # These need to be liftover before 
        replacement_snp_ids = convert_chr_positions_to_rsids(to_fix$SNP)
        single_eqtl[str_detect(single_eqtl$SNP, ':')]$SNP=replacement_snp_ids$rsid
    }

    

    #eQTL col rename
    # single_eqtl <- dplyr::rename(single_eqtl,
    #   p_value = p,
    #   variant_id = SNP,
    #   standard_error = se
    # )
    return (single_eqtl)
}

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
    # 'P' = "p_value",
    # 'P-value' = "p_value",
    'P' = "p",
    'P-value' = "p",
    'REF' = "reference_allele",
    'ALT' = "alternative_allele",

    'Allele1' = "effect_allele",
    'A1' = "effect_allele",
    'A_EFF' = "effect_allele",
    'EA' = "effect_allele",

    'Allele2' = "other_allele",
    'AX' = 'other_allele',
    'A_NONEFF' = 'other_allele',
    'NEA' = 'other_allele',

    'CHROM' = "chromosome",
    'CHR' = "chromosome",

    # 'ID' = "variant_id",
    # 'MarkerName' = "variant_id",
    # 'SNP' = "variant_id",
    # 'rsID' = "variant_id",
    'ID' = "SNP",
    'MarkerName' = "SNP",
    'rsID' = "SNP",
    'POS' = "base_pair_location",
    'pos' = "base_pair_location",
    'BP' = "base_pair_location",

    'BETA' = "beta",
    'Effect' = "beta",

    # 'SE' = "standard_error",
    # 'StdErr' = "standard_error",
    'SE' = "se",
    'StdErr' = "se",
    'Freq1' = 'eaf',
    'A1_FREQ' = 'eaf',
    'Freq_EFF' = 'eaf',
    'EAF' = 'eaf',

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

convert_chr_positions_to_rsids  <- function(variant_positions){
  # Since most of the GWAS summary statistics come with an rsid denotion
  # we should make sure that we mapp these ids to the rsids.
  # variant_positions = to_fix$SNP
  Data2 = data.frame(str_split_fixed(variant_positions, '[:_]', 4))
  Data2$X1 <-str_replace(Data2$X1, "chr", "")
  Data2$X3_2=as.numeric(Data2$X2)-1
  Data2$original_id=variant_positions
  Data2_t=Data2[c('X1','X3_2','X2','original_id')]
  Data2_t$R = paste0(Data2_t$X1,':',Data2_t$X2)

  # before performing this we may want to liftover any hg37 positions to hg38
  
  write.table(Data2_t, file=paste0("tmp_map.bed"), sep = "\t", quote = FALSE, row.names = FALSE,col.names=FALSE)
  bin = "bcftools query -f'%CHROM-%POS\t%ID\t%REF\t%ALT\n' -R tmp_map.bed rsid_vcf_with_ref.vcf.gz > mappings.tsv"
  # logfile <- tempfile()
  
  rc <- system(bin)
  Data2_t$V1=paste0(Data2_t$X1,'-',Data2_t$X2)
  rownames(Data2_t)=Data2_t$map1
  Data2_t$rsid_ref_alt = paste0(Data2_t$V1,'-',Data2$X3,'-',Data2$X4)
  Data2_t$rsid_alt_ref = paste0(Data2_t$V1,'-',Data2$X4,'-',Data2$X3)
  mappings=fread('mappings.tsv',sep="\t",header=FALSE)
  mappings <- mappings %>% separate_rows(V4)
  mappings <- mappings %>% separate_rows(V3)
  mappings$original_id <- 'original_id'
  mappings$rsid_ref_alt = paste0(mappings$V1,'-',mappings$V3,'-',mappings$V4)
  mappings$rsid_alt_ref = paste0(mappings$V1,'-',mappings$V4,'-',mappings$V3)
  Data2_t$rsid = Data2_t$original_id

  merged2 = left_join(Data2_t, mappings, by = join_by(rsid_ref_alt == rsid_ref_alt),  multiple = "any",unmatched = "drop")
  merged2 = merged2[!is.na(merged2$V2),]
  merged2 = subset(merged2, select = -c(rsid) )
  names(merged2)[names(merged2) == "V2"] <- "rsid"
  merged3 = merged2[c('rsid','rsid_ref_alt')]  

  merged32 = left_join(Data2_t, mappings, by = join_by(rsid_ref_alt == rsid_alt_ref),  multiple = "any",unmatched = "drop")
  merged32 = merged32[!is.na(merged32$V2),]
  merged32 = subset(merged32, select = -c(rsid) )
  names(merged32)[names(merged32) == "V2"] <- "rsid"
  merged32 = merged32[c('rsid','rsid_alt_ref')]   
  # names(merged32)[names(merged32) == "rsid_alt_ref"] <- "rsid_ref_alt" 

  Data2_t = Data2_t %>% 
    rows_update(unique(merged32), by = "rsid_alt_ref") 

  Data2_t = Data2_t %>% 
    rows_update(unique(merged3), by = "rsid_ref_alt")

  file.remove("mappings.tsv")
  file.remove("tmp_map.bed")
  return(Data2_t)
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

add_position <- function (df, variant_colname = 'SNP', genotype_prefix){
    stopifnot(variant_colname %in% colnames(df))
    bim_filename <- paste(genotype_prefix, 'bim', sep='.')
    bim <- fread(bim_filename, col.names = c('chromosome', 'rsid', 'base_pair_location'), select = c(1, 2, 4))
    m <- dplyr::inner_join(df, bim, by = setNames('rsid', variant_colname))
    return(m)
}
