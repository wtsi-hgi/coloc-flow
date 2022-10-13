
make_cojo_df <- function(GWAS_matched_SNPS_with_eQTL){
          GWAS_matched_SNPS_with_eQTL$FREQ=freqs$FreqA1[match(GWAS_matched_SNPS_with_eQTL$variant_id,freqs$SNP)]
          idx=which(GWAS_matched_SNPS_with_eQTL$effect_allele!=freqs$effect_allele)
          GWAS_matched_SNPS_with_eQTL$FREQ[idx]=1-GWAS_matched_SNPS_with_eQTL$FREQ[idx]
          # GWAS_matched_SNPS_with_eQTL$se=sqrt(GWAS_matched_SNPS_with_eQTL$varbeta)
          Cojo_Dataframe=GWAS_matched_SNPS_with_eQTL[,c("variant_id","effect_allele","other_allele","FREQ","beta","standard_error","p_value","N")]
          names(Cojo_Dataframe)=c("SNP" , "A1" ,  "A2"  , "freq", "b"  ,  "se" ,  "p" ,   "N")
          # from https://yanglab.westlake.edu.cn/software/gcta/#COJO
          # A1 -- the effect allele
          # A2 -- the other allele
          # freq -- frequency of the effect allele
          # b -- effect size
          # p -- p-value
          # N -- sample size
          Cojo_Dataframe=Cojo_Dataframe[!is.na(Cojo_Dataframe$freq)]# <- 0
          Cojo_Dataframe=Cojo_Dataframe[!Cojo_Dataframe$N <10]# <- 10 #Cojo doesnt like sample sizes smaller than 10
          Cojo_Dataframe = transform(Cojo_Dataframe, freq = as.numeric(freq), N = as.numeric(N),b = as.numeric(b))
    return (Cojo_Dataframe)
}

# call gcta program
run_gcta <- function (bin='gcta', args){
    print(paste(bin, args, collapse = ' '))
    rc <- system2(command = bin, args = args)
    stopifnot(rc == 0)
}

# call gcta program in COJO-mode
run_cojo <- function (bin=NULL, bfile, marker_list, summary_stat, out_prefix){
    gcta_args <- c(
        '--bfile', bfile,
        '--extract', marker_list,     # limit the COJO analysis in a certain genomic region
        '--cojo-slct',                # select independently associated SNPs
        '--cojo-p', '1e-4',           # p-value to declare a genome-wide significant hit
        '--cojo-file', summary_stat,  # summary-level statistics from a GWAS
        '--out', out_prefix
    )

    cmd <- list(args = gcta_args)
    if(!is.null(bin)){
        cmd$bin <- bin
    }
    do.call(run_gcta, cmd)

    outfile <- paste0(out_prefix, ".jma.cojo")
    return(outfile)
}
