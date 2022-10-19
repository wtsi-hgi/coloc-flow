# https://yanglab.westlake.edu.cn/software/gcta/#COJO
library(data.table)
requireNamespace(dplyr)

make_cojo_df <- function(GWAS_matched_SNPS_with_eQTL, freqs){
    GWAS_matched_SNPS_with_eQTL$FREQ = freqs$FreqA1[match(GWAS_matched_SNPS_with_eQTL$variant_id,freqs$SNP)]
    idx = which(GWAS_matched_SNPS_with_eQTL$effect_allele!=freqs$effect_allele)
    GWAS_matched_SNPS_with_eQTL$FREQ[idx] = 1-GWAS_matched_SNPS_with_eQTL$FREQ[idx]
    # GWAS_matched_SNPS_with_eQTL$se=sqrt(GWAS_matched_SNPS_with_eQTL$varbeta)

    Cojo_Dataframe <- dplyr::select(GWAS_matched_SNPS_with_eQTL,
        SNP=variant_id,
        A1=effect_allele,  # the effect allele
        A2=other_allele,   # the other allele
        freq=FREQ,         # frequency of the effect allele
        b=beta,            # effect size
        se=standard_error,
        p=p_value,         # p-value
        N                  # sample size
    )
    Cojo_Dataframe = Cojo_Dataframe[!is.na(Cojo_Dataframe$freq)]  # <- 0
    Cojo_Dataframe = Cojo_Dataframe[!Cojo_Dataframe$N < 10]  # Cojo doesnt like sample sizes smaller than 10
    Cojo_Dataframe = transform(Cojo_Dataframe, freq = as.numeric(freq), N = as.numeric(N), b = as.numeric(b))
    return (Cojo_Dataframe)
}

# a small wrapper for system call with return code check
run_tool <- function (bin, args){
    print(paste(c(bin, args), collapse = ' '))
    rc <- system2(command = bin, args = args)
    stopifnot(rc == 0)
}

# call gcta program
run_gcta <- function (bin = NULL, args){
    if (is.null(bin)) bin <- 'gcta'
    run_tool(bin = bin, args = args)
}

# call gcta program in COJO-mode
# bfile -- path plink file with LD reference panel
run_cojo <- function (bin = NULL, bfile, marker_list = NULL, conditional_markers = NULL, summary_stat, out_prefix){
    gcta_args <- c(
        '--bfile', bfile,
        '--cojo-slct',                # select independently associated SNPs
        '--cojo-p', '1e-4',           # p-value to declare a genome-wide significant hit
        '--cojo-file', summary_stat,  # summary-level statistics from a GWAS
        '--out', out_prefix
    )

    if (!is.null(conditional_markers)){
        gcta_args <- append(gcta_args,
            c('--cojo-cond', conditional_markers)  # analysis conditional on the given list of SNPs
        )
    }

    if (!is.null(marker_list)){
        gcta_args <- append(gcta_args,
            c('--extract', marker_list)  # limit the COJO analysis in a certain genomic region
        )
    }

    run_gcta(bin = bin, args = gcta_args)

    out <- list(
        independent_signals = paste0(out_prefix, ".jma.cojo"),
        conditional_analysis = paste0(out_prefix, ".cma.cojo")
    )
    return(out)
}

run_cojo_on_locus <- function (gcta_bin = NULL, plink2_bin = NULL,
                               bfile, chrom, start, end,
                               conditional_markers = NULL, summary_stat, freqs, out_prefix){
    variants_of_interest <- dplyr::filter(summary_stat,
        chromosome == chrom,
        base_pair_location >= start,
        base_pair_location <= end
    )

    Cojo_Dataframe <- make_cojo_df(variants_of_interest, freqs = freqs)

    cojo_filename <- paste0(variant_id, '_', GWAS_name, "_sum.txt")
    fwrite(Cojo_Dataframe, file = cojo_filename, row.names = F, quote = F, sep = "\t")

    locus_filename <- paste(basename(bfile), chrom, start, end, sep = '-')
    extract_locus(
        bin = plink2_bin,
        genotypes_prefix = bfile,
        chrom = chrom, start = start, end = end,
        out_prefix = locus_filename
    )

    cojo <- run_cojo(
        bin = gcta_bin,
        bfile = locus_filename,
        marker_list = marker_filename,
        summary_stat = cojo_filename,
        out_prefix = out_prefix
    )
    return(cojo)
}

# call plink program
run_plink2 <- function (bin = NULL, args){
    if (is.null(bin)) bin <- 'plink2'
    run_tool(bin = bin, args = args)
}

# extracts genomic interval from bfile using plink2
extract_locus <- function (bin = NULL, genotypes_prefix, chrom, start, end, out_prefix){
    range_filename <- tempfile(tmpdir = '.', fileext = '.txt')
    df <- data.table(chrom = chrom, start = start, end = end, label = 'locus')
    fwrite(df, range_filename, sep = '\t', col.names = F)

    plink_args <- c(
        '--bfile', genotypes_prefix,
        '--extract', 'bed1', range_filename,
        '--make-bed',
        '--out', out_prefix
    )

    run_plink2(bin = bin, args = plink_args)
    file.remove(range_filename)

    return(out_prefix)
}
