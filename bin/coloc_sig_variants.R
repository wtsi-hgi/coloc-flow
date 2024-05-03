#!/usr/bin/env Rscript
library(data.table)
library(coloc)
# library(susieR)
library(optparse)
library(magrittr)
library(dplyr)
library(tidyr)
requireNamespace('dplyr')
requireNamespace('tidyr')
library("stringr") 
eqtl_significance_threshold <- 5e-5
gwas_significance_threshold <- 5e-7
# --config /lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline_ip13/input.yml
option_list <- list(
    make_option('--gwas', action="store", help="path to GWAS summary statistic"),
    make_option('--rs', action="store", help="variant ID"),
    make_option('--bfile', action="store", help="path to plink genotypes prefix"),
    make_option('--eqtl', action="store", help="path to eqtl data"),
    make_option('--eqtl_snps', action="store", help = "path to eqtl snp_pos.txt file"),
    make_option('--eqtl_samples', action="store", default=192, help = "number of samples in eQTL study"),
    make_option('--config', action="store", defaul=NULL, help="path to yaml-file with data configs"),
    make_option('--plink2_bin', action="store", default=NULL, help = "path to custom plink2 binary"),
    make_option('--gcta_bin', action="store", default=NULL, help = "path to custom gcta binary")
)
args <- parse_args(OptionParser(option_list=option_list))

source('dataIO.R')
source('cojo.R')
source('helpers.R')

eQTL = args$eqtl
eqtl_marker_file = args$eqtl_snps
eqtl_samples_number = args$eqtl_samples
GWAS = args$gwas
variant = args$rs
bfile = args$bfile
plink2_bin = args$plink2_bin
gcta_bin = args$gcta_bin
contig = args$config


# ---gwas GIGASTROKE_CES_EUR_hg19_harmonised.tsv.gz             --rs rs3756011             --bfile plink_genotypes             --eqtl Astrocytes.4.gz             --eqtl_snps snp_pos.txt             --config /lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline_ip13/input.yml
# eQTL = 'Astrocytes.4.gz'
# eqtl_marker_file = 'snp_pos.txt'
# eqtl_samples_number = 192
# GWAS = 'GIGASTROKE_CES_EUR_hg19_harmonised.tsv.gz'
# variant = 'rs3756011'
# bfile = './plink_genotypes'
# plink2_bin = NULL
# gcta_bin = NULL
# contig = '/lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline_ip13/input.yml'
# args$eqtl_snps = 'snp_pos.txt'

GWAS_name = tools::file_path_sans_ext(basename(GWAS), compression = T)
eQTL_name = strsplit(tools::file_path_sans_ext(basename(eQTL)), "\\.")[[1]]
eQTL_name = paste(eQTL_name[1:length(eQTL_name)-1], collapse = "_")

config <- read_config(contig, eQTL_name)
if (is.null(config$build)){
    bd ='hg38'
}else{
    bd =config$build
}

config <- read_config(contig, GWAS_name)
if (is.null(config$build)){
    bd2 ='hg38'
}else{
    bd2 =config$build
}

eqtl_marker_data <- read_eqtl_marker_file(args$eqtl_snps, build = bd2) # Dependant on which build we are using we preload the SNP variant file to ease the mapping between different builds.
Full_GWAS_Sum_Stats = load_eqtl(GWAS, marker.data = eqtl_marker_data, build = bd2)
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats,threshold=gwas_significance_threshold)

eqtl_marker_data <- read_eqtl_marker_file(args$eqtl_snps, build = bd) 
single_eqtl1 = load_eqtl(eQTL, marker.data = eqtl_marker_data, build = bd)

row1 = Significant_GWAS_Signals[Significant_GWAS_Signals$SNP == variant,]
base_pair_location <- row1[["base_pair_location"]]
SNP <- row1[["SNP"]]
chromosome1 <- row1[["chromosome"]]
print(paste('Running GWAS variant', SNP))

locus_start <- base_pair_location - 1e6
locus_end <- base_pair_location + 1e6

variants_of_interest <- dplyr::filter(Full_GWAS_Sum_Stats,
    chromosome == chromosome1,
    between(base_pair_location, locus_start, locus_end)
)

gwas_build = bd2
# if(!is.null(config$build)){
#     gwas_build <- config$build
# } else{
#     gwas_build <- get_df_build_version(df = variants_of_interest)
# }

# if(gwas_build != 'hg38'){
#     message(paste('Convert GWAS positions from', gwas_build, 'to hg38'))
#     ch <- load_chain_file(from = gwas_build, to = 'hg38')
#     variants_of_interest <- lift_over_df(variants_of_interest, chain = ch)
#     variants_of_interest <- dplyr::filter(variants_of_interest, chromosome == chromosome1)

#     locus_start <- min(variants_of_interest$base_pair_location)
#     locus_end <- max(variants_of_interest$base_pair_location)
# }

# Here we have to add the rsids if these are not present in the gwas file. rs28887923

locus_filename <- paste(basename(bfile), chromosome1, locus_start, locus_end, sep = '-')
extract_locus(
    bin = plink2_bin,
    genotypes_prefix = bfile,
    chrom = chromosome1,
    start = locus_start,
    end = locus_end,
    filters = c('--maf', '0.01'),
    out_prefix = locus_filename
)

if(!'eaf' %in% colnames(variants_of_interest)){
    if(!is.null(config$frequency_not_available) && config$frequency_not_available){
        message('Calculating allele frequencies from LD panel')
        variants_of_interest %<>% add_freq(
            freq = get_plink_freq(plink2_bin, locus_filename)
        )
    } else {
        text <- paste('No AF column in GWAS. If you wish to use AF of LD-panel',
                      'please use `frequency_not_available` in yaml-config')
        stop(text)
    }
}

if(!is.null(config$type)){
    gwas_type <- config$type
    stopifnot(gwas_type %in% c('cc', 'quant'))
} else {
    stop('Please provide GWAS type (cc or quant) in yaml-config')
}

if(!'N' %in% colnames(variants_of_interest)){
    if('samples_number' %in% names(config)){
        variants_of_interest %<>% dplyr::mutate(N = config$samples_number)
    } else {
        stop('No information about number of samples in GWAS. Please provide it in yaml-config')
    }
}

Cojo_Dataframe <- make_cojo_df(variants_of_interest)

genes_of_interest <- dplyr::filter(single_eqtl1,
    chromosome == chromosome1,
    between(base_pair_location, locus_start, locus_end)
)
eqtl_genes <- unique(genes_of_interest$gene)

cojo_filename <- paste0(SNP, '_', GWAS_name, "_sum.txt")
fwrite(Cojo_Dataframe, file = cojo_filename, row.names = F, quote = F, sep = "\t")

cojo_out <- run_cojo(
    bin = gcta_bin,
    bfile = locus_filename,
    summary_stat = cojo_filename,
    pvalue = cojo_strict_threshold,
    out_prefix = paste(SNP, GWAS_name, "step1", sep = "_")
)
# bin = NULL, bfile, marker_list = NULL, conditional_markers = NULL, summary_stat, pvalue, out_prefix
coloc_results <- list()
if(!file.exists(cojo_out$independent_signals)){
    warning('No independant SNPs in GWAS found')
    independent_SNPs <- dplyr::rename(row1, SNP = SNP, p = p)
} else {
    independent_SNPs <- fread(cojo_out$independent_signals)
}
independent_signals <- dplyr::filter(independent_SNPs, p < gwas_signal_threshold)$SNP

for( GWAS_signal in independent_signals){
    all_but_one <- setdiff(independent_SNPs$SNP, GWAS_signal)

    if(length(all_but_one) > 0){
        independent_markerfile <- paste0(GWAS_signal, "_independent.snp")
        writeLines(all_but_one, con = independent_markerfile)

        # some SNPs have large difference of allele frequency between the GWAS summary data and the reference sample, hence are removed.
        cojo_cond_out <- run_cojo(
            bin = gcta_bin,
            bfile = locus_filename,
            summary_stat = cojo_filename,
            pvalue = cojo_strict_threshold,
            conditional_markers = independent_markerfile,
            out_prefix = paste(SNP, GWAS_name, "step2", sep = '_')
        )

        conditioned_dataset <- fread(cojo_cond_out$conditional_analysis)
        gwas_data <- conditioned_dataset
    } else {
        gwas_data <- variants_of_interest
    }

    for (qtl1 in eqtl_genes){

        single_eqtl = genes_of_interest[genes_of_interest$gene==qtl1, ]
    
        single_eqtl2 = single_eqtl
        rownames(single_eqtl2) <- single_eqtl2$SNP
        single_eqtl2$N <- eqtl_samples_number

        if (min(single_eqtl2$p) < eqtl_significance_threshold){
            print('yes')
        
            Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2, source = 'eqtl')

    
            eqtl_summary_file <- paste(SNP, qtl1, eQTL_name, "eqtl_sum.txt", sep = "_")
            fwrite(Cojo_Dataframe_eqtl, file = eqtl_summary_file, row.names = F, quote = F, sep = "\t")

            eqtl_cojo_out <- run_cojo(
                bin = gcta_bin,
                bfile = locus_filename,
                summary_stat = eqtl_summary_file,
                pvalue = cojo_strict_threshold,
                out_prefix = paste(SNP, qtl1, eQTL_name, "eqtl_step1", sep = "_")
            )
            if(!file.exists(eqtl_cojo_out$independent_signals)){
                warning('No independant SNPs in eQTL found')
                independent_SNPs_eQTL <- dplyr::rename(single_eqtl2, SNP = SNP, p = p) %>%
                  dplyr::slice_min(p, n = 1, with_ties = F)
            } else {
                independent_SNPs_eQTL <- fread(eqtl_cojo_out$independent_signals)
            }

            for( independent_eqtl_SNP_to_contition_on in independent_SNPs_eQTL$SNP){
                # print(independent_eqtl_SNP_to_contition_on)}
                all_but_one_eqtl <- setdiff(independent_SNPs_eQTL$SNP, independent_eqtl_SNP_to_contition_on)
            
                if(length(all_but_one_eqtl) > 0){
                    eqtl_independant_markerfile <- paste0(independent_eqtl_SNP_to_contition_on, "_eqtl_independent.snp")
                    writeLines(all_but_one_eqtl, con = eqtl_independant_markerfile)

                    eqtl_cond_cojo_out <- run_cojo(
                        bin = gcta_bin,
                        bfile = locus_filename,
                        summary_stat = eqtl_summary_file,
                        pvalue = cojo_strict_threshold,
                        conditional_markers = eqtl_independant_markerfile,
                        out_prefix = paste(SNP, qtl1, eQTL_name, independent_eqtl_SNP_to_contition_on, "eqtl_step2", sep = "_")
                    )

                    conditioned_dataset_eQTL <- fread(eqtl_cond_cojo_out$conditional_analysis)
                    eqtl_data <- conditioned_dataset_eQTL
                } else {
                    eqtl_data <- single_eqtl2
                }

                # eqtl_data <- align_datasets(gwas_data, eqtl_data)
                D1 <- prepare_coloc_table(gwas_data)
                D2 <- prepare_coloc_table(eqtl_data)

                D1_l <- prepare_coloc_list(D1, N = first(Cojo_Dataframe$N), type = gwas_type)
                D2_l <- prepare_coloc_list(D2, N = eqtl_samples_number)

                # FIXME align effect alleles in both datasets
                colo.res <- coloc.abf(D1_l, D2_l)
                colo_res=data.frame(t(colo.res$summary))
                print('|||Coloc result:|||')
                print(colo_res)

                if (colo_res$PP.H4.abf > coloc_threshold){

                    fig1_filename <- paste(GWAS_name, GWAS_signal, eQTL_name, qtl1, independent_eqtl_SNP_to_contition_on, chromosome1, locus_start, locus_end, 'coloc.jpg', sep='_')
                    fig2_filename <- paste(GWAS_name, GWAS_signal, eQTL_name, qtl1, independent_eqtl_SNP_to_contition_on, chromosome1, locus_start, locus_end, 'rplot.jpg', sep='_')

                    colo_df <- data.table(
                      gwas_name = GWAS_name,
                      gwas_input = variant,
                      locus = paste(chromosome1, locus_start, locus_end, sep = '-'),
                      gwas_lead = GWAS_signal,
                      eqtl_name= eQTL_name,
                      gene = qtl1,
                      eqtl_lead = independent_eqtl_SNP_to_contition_on,
                      pp_h4 = colo.res$summary[['PP.H4.abf']],
                      colocolised_snp = subset(colo.res$results, SNP.PP.H4>0.01)$snp,
                      figure_data = file.path(getwd(), fig2_filename),
                      figure_coloc = file.path(getwd(), fig1_filename)
                    )

                    coloc_results <- append(coloc_results, list(colo_df))

                    try({
                        jpeg(fig1_filename)
                            sensitivity(colo.res, paste0("H4 > ", coloc_threshold))
                        dev.off()
                    }, silent = F)

                    if(nrow(colo_res)>0){
                        p1 <- plot_ggwas(D1, position, pvalues, xlim=c(locus_start, locus_end),
                                         snp_column = 'snp', highlight_snps = colo_df$colocolised_snp) +
                            ggplot2::labs(title = GWAS_name, subtitle = GWAS_signal)
                        p2 <- plot_ggwas(D2, position, pvalues, xlim=c(locus_start, locus_end),
                                         snp_column = 'snp', highlight_snps = colo_df$colocolised_snp) +
                            ggplot2::labs(title = paste(eQTL_name, qtl1, sep = ' · '), subtitle = independent_eqtl_SNP_to_contition_on)
                        p <- gridExtra::grid.arrange(p1, p2, nrow = 2, ncol = 1,
                                                     top = grid::textGrob(paste('chromosome', chromosome1)))
                        ggplot2::ggsave(plot = p, filename = fig2_filename)
                    }
                }
            }
        }
    }
}

if(length(coloc_results) > 0){
    coloc_df <- rbindlist(coloc_results)
    fwrite(coloc_df, 'coloc_results.csv')
}
