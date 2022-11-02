requireNamespace('biomaRt')

gwas_significance_threshold <- 5e-8
eqtl_significance_threshold <- 5e-5

get_gwas_significant_signals <- function (df, threshold=gwas_significance_threshold){
    df_sign <- df[df$p_value < threshold]
    return(df_sign)
}


ensembl38 <- biomaRt::useEnsembl(biomart = "snps", dataset = 'hsapiens_snp')
ensembl37 <- biomaRt::useEnsembl(biomart = "snps", dataset = 'hsapiens_snp', version = 'GRCh37')

# retrieve SNP position from BioMart
get_snp_position <- function (rs, mart){
    filter_values <- c('snp_filter', 'snp_synonym_filter')
    for(filter_value in filter_values){
        response <- biomaRt::getBM(
            mart = mart,
            attributes = c('chr_name', 'chrom_start', 'allele'),
            filters = filter_value,
            values = rs
        )
        if(nrow(response != 0)){
            break
        }
    }
    if(nrow(response) > 1){
        response <- dplyr::filter(response, chr_name %in% c(1:22, 'X', 'Y', 'XY', 'M', 'MT'))
    }
    position <- as.integer(response$chrom_start)
    alleles <- unlist(strsplit(response$allele, '/'))
    if(any(nchar(alleles) > 1)){
        position <- position - 1  # in VCF indels have position one nucleotide less
    }
    return(position)
}

# reveals SNP position build version (hg37/hg38)
get_snp_build_version <- function (rs, pos, mart37 = ensembl37, mart38 = ensembl38){
    pos37 <- get_snp_position(rs, mart = mart37)
    pos38 <- get_snp_position(rs, mart = mart38)
    stopifnot(pos37 != pos38)

    if(pos == pos37)
        return('hg19')

    if(pos == pos38)
        return('hg38')

    return(NULL)
}

# reveals dataframe build version (hg37/hg38)
get_df_build_version <- function (df, mart37 = ensembl37, mart38 = ensembl38){
    df_rs <- dplyr::filter(df,
      startsWith(variant_id, "rs"),
      nchar(effect_allele) == 1,
      nchar(other_allele) == 1
    )
    stopifnot(nrow(df_rs) > 0)
    row <- dplyr::slice_sample(df_rs, n = 1)
    message(paste('Chosing', row$variant_id, 'to reveal genome build'))
    get_snp_build_version(rs = row$variant_id, pos = row$base_pair_location,
                          mart37 = mart37, mart38 = mart38)
}

# downloads and loads liftOver chain-file
load_chain_file <- function (from = 'hg19', to = 'hg38'){
    filename <- paste0(from, "To", stringr::str_to_title(to), ".over.chain")
    if(!file.exists(filename)){
        gz_filename <- paste(filename, 'gz', sep = '.')
        message(paste("Downloading", gz_filename))
        url <- file.path("https://hgdownload.cse.ucsc.edu/goldenpath", from, "liftOver", gz_filename)
        curl::curl_fetch_disk(url, gz_filename)
        R.utils::gunzip(gz_filename, remove = T)
    }
    rtracklayer::import.chain(filename)
}

lift_over_bp <- function (chrom, bp, chain){
    irange <- IRanges::IRanges(start = bp, end = bp)
    grange <- GenomicRanges::GRanges(seqnames = chrom, ranges = irange)
    GenomeInfoDb::seqlevelsStyle(grange) <- "UCSC"

    grange_new <- unlist(rtracklayer::liftOver(grange, chain))
    bp_new <- BiocGenerics::start(grange_new)
    return(bp_new)
}

lift_over_df <- function (df, chain){
    stopifnot(all(c('chromosome', 'base_pair_location') %in% colnames(df)))
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df,
      keep.extra.columns = T,
      seqnames.field = 'chromosome',
      start.field = 'base_pair_location',
      end.field = 'base_pair_location'
    )
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

    gr_new <- unlist(rtracklayer::liftOver(gr, chain))
    df_new <- data.table::as.data.table(gr_new)
    df_new$seqnames <- as.integer(gsub('chr', '', df_new$seqnames))
    df_new <- dplyr::rename(df_new, base_pair_location = start, chromosome = seqnames)
    df_new <- dplyr::select(df_new, -end, -width, -strand)
    return(df_new)
}

thin_out_gwas_data <- function(df, threshold = 1e-5){
    sign_snps <- dplyr::filter(df, p_value <= threshold)
    other_snps <- dplyr::filter(df, p_value > threshold)
    other_snps <- dplyr::slice_sample(other_snps, prop=0.1)
    rbind(sign_snps, other_snps)
}

plot_gwas <- function (df){
    gwas <- thin_out_gwas_data(df)
    qqman::manhattan(gwas, chr = 'chromosome', bp = 'base_pair_location', snp = 'variant_id', p = 'p_value')
}
