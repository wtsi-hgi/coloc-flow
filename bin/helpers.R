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
    position <- biomaRt::getBM(
        mart = mart,
        attributes = 'chrom_start',
        filters = 'snp_filter',
        values = rs
    )
    return(as.integer(position))
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
