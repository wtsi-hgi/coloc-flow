gwas_significance_threshold <- 5e-8
eqtl_significance_threshold <- 5e-5

get_gwas_significant_signals <- function (df, threshold=gwas_significance_threshold){
    df_sign <- df[df$p_value < threshold]
    return(df_sign)
}
