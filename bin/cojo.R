# https://yanglab.westlake.edu.cn/software/gcta/#COJO
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

# call gcta program
run_gcta <- function (bin='gcta', args){
    print(paste(c(bin, args), collapse = ' '))
    rc <- system2(command = bin, args = args)
    stopifnot(rc == 0)
}

# call gcta program in COJO-mode
run_cojo <- function (bin = NULL, bfile, marker_list, conditional_markers = NULL, summary_stat, out_prefix){
    gcta_args <- c(
        '--bfile', bfile,
        '--extract', marker_list,     # limit the COJO analysis in a certain genomic region
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

    cmd <- list(args = gcta_args)
    if(!is.null(bin)){
        cmd$bin <- bin
    }
    do.call(run_gcta, cmd)

    out <- list(
        independent_signals = paste0(out_prefix, ".jma.cojo"),
        conditional_analysis = paste0(out_prefix, ".cma.cojo")
    )
    return(out)
}
