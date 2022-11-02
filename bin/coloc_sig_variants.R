#!/usr/bin/env Rscript
library(data.table)
library(coloc)
# library(susieR)
library(optparse)

option_list <- list(
    make_option('--gwas', action="store", help="path to GWAS summary statistic"),
    make_option('--rs', action="store", help="variant ID"),
    make_option('--bfile', action="store", help="path to plink genotypes prefix"),
    make_option('--eqtl', action="store", help="path to eqtl data"),
    make_option('--eqtl_snps', action="store", help = "path to eqtl snp_pos.txt file"),
    make_option('--eqtl_samples', action="store", default=192, help = "number of samples in eQTL study"),
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

GWAS_name = tools::file_path_sans_ext(basename(GWAS))
eQTL_name = strsplit(tools::file_path_sans_ext(basename(eQTL)), "\\.")[[1]]
eQTL_name = paste(eQTL_name[1:length(eQTL_name)-1], collapse = "_")

# freq_file=paste0(GWAS_name, ".frqx")
# freqs = read_freqs(freq_file)

Full_GWAS_Sum_Stats = load_GWAS(GWAS)$map
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats)
single_eqtl1 = load_eqtl(eQTL, eqtl_marker_file)

row1 = Significant_GWAS_Signals[Significant_GWAS_Signals$variant_id == variant,]
base_pair_location <- row1[["base_pair_location"]]
variant_id <- row1[["variant_id"]]
chromosome1 <- row1[["chromosome"]]
print(paste('Running GWAS variant', variant_id))

locus_start <- base_pair_location - 1e6
locus_end <- base_pair_location + 1e6

variants_of_interest <- dplyr::filter(Full_GWAS_Sum_Stats,
    chromosome == chromosome1,
    between(base_pair_location, locus_start, locus_end)
)

variant_build <- get_df_build_version(df = variants_of_interest)
if(variant_build != 'hg38'){
    message(paste('Convert GWAS positions from', variant_build, 'to hg38'))
    ch <- load_chain_file(from = variant_build, to = 'hg38')
    variants_of_interest <- lift_over_df(variants_of_interest, chain = ch)

    locus_start <- min(variants_of_interest$base_pair_location)
    locus_end <- max(variants_of_interest$base_pair_location)
}

Cojo_Dataframe <- make_cojo_df(variants_of_interest)

genes_of_interest <- dplyr::filter(single_eqtl1,
    chromosome == chromosome1,
    between(base_pair_location, locus_start, locus_end)
)
all_unique_eQTL_signals_in_this_GWAS_range = unique(genes_of_interest$gene)

cojo_filename <- paste0(variant_id, '_', GWAS_name, "_sum.txt")
fwrite(Cojo_Dataframe, file = cojo_filename, row.names = F, quote = F, sep = "\t")

locus_filename <- paste(basename(bfile), chromosome1, locus_start, locus_end, sep = '-')
extract_locus(
    bin = plink2_bin,
    genotypes_prefix = bfile,
    chrom = chromosome1,
    start = locus_start,
    end = locus_end,
    out_prefix = locus_filename
)

cojo_out <- run_cojo(
    bin = gcta_bin,
    bfile = locus_filename,
    summary_stat = cojo_filename,
    pvalue = gwas_significance_threshold,
    out_prefix = paste(variant_id, GWAS_name, "step1", sep = "_")
)

independent_SNPs = fread(cojo_out$independent_signals)

coloc_results <- list()
for( i_GWAS in 1:nrow(independent_SNPs)){
    GWAS_signal = independent_SNPs[i_GWAS]$SNP
    independent_markerfile <- paste0(GWAS_signal, "_independent.snp")
    writeLines(GWAS_signal, con = independent_markerfile)
    cojo_cond_out <- run_cojo(
        bin = gcta_bin,
        bfile = locus_filename,
        summary_stat = cojo_filename,
        pvalue = gwas_significance_threshold,
        conditional_markers = independent_markerfile,
        out_prefix = paste(variant_id, GWAS_name, "step2", sep = '_')
    )
    conditioned_dataset <- combine_cojo_results(independent_signals = cojo_out$independent_signals,
                                                conditional_signals = cojo_cond_out$conditional_analysis,
                                                lead_snp = GWAS_signal)
    # some SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample, hence are removed.
    D1 <- prepare_coloc_table(conditioned_dataset)

    for (qtl1 in all_unique_eQTL_signals_in_this_GWAS_range){
        single_eqtl = genes_of_interest[genes_of_interest$gene==qtl1, ]
        single_eqtl2 = single_eqtl
        rownames(single_eqtl2) <- single_eqtl2$SNP
        single_eqtl2$N <- eqtl_samples_number

        if (min(single_eqtl2$p) < eqtl_significance_threshold){
            Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2, source = 'eqtl')

            eqtl_summary_file <- paste(variant_id, qtl1, eQTL_name, "eqtl_sum.txt", sep = "_")
            fwrite(Cojo_Dataframe_eqtl, file = eqtl_summary_file, row.names = F, quote = F, sep = "\t")

            eqtl_cojo_out <- run_cojo(
                bin = gcta_bin,
                bfile = locus_filename,
                summary_stat = eqtl_summary_file,
                pvalue = eqtl_significance_threshold,
                out_prefix = paste(variant_id, qtl1, eQTL_name, "eqtl_step1", sep = "_")
            )
            independent_SNPs_eQTL <- fread(eqtl_cojo_out$independent_signals)

            for( i_eQTL in 1:nrow(independent_SNPs_eQTL)){
                print(i_eQTL)
                independent_eqtl_SNP_to_contition_on = independent_SNPs_eQTL[i_eQTL]$SNP
                eqtl_independant_markerfile <- paste0(independent_eqtl_SNP_to_contition_on, "_eqtl_independent.snp")
                writeLines(independent_eqtl_SNP_to_contition_on, con = eqtl_independant_markerfile)

                eqtl_cond_cojo_out <- run_cojo(
                    bin = gcta_bin,
                    bfile = locus_filename,
                    summary_stat = eqtl_summary_file,
                    pvalue = eqtl_significance_threshold,
                    conditional_markers = eqtl_independant_markerfile,
                    out_prefix = paste(variant_id, qtl1, eQTL_name, independent_eqtl_SNP_to_contition_on, "eqtl_step2", sep = "_")
                )

                conditioned_dataset_eQTL <- combine_cojo_results(
                    independent_signals = eqtl_cojo_out$independent_signals,
                    conditional_signals = eqtl_cond_cojo_out$conditional_analysis,
                    lead_snp = independent_eqtl_SNP_to_contition_on
                )

                D2 <- prepare_coloc_table(conditioned_dataset_eQTL)

                SNPs_To_Colocalise = intersect(D1$snp, D2$snp)
                D1_col = D1[D1$snp %in% SNPs_To_Colocalise, ]
                D2_col = D2[D2$snp %in% SNPs_To_Colocalise, ]

                D1_l <- prepare_coloc_list(D1_col, N = first(Cojo_Dataframe$N))
                D2_l <- prepare_coloc_list(D2_col, N = eqtl_samples_number)

                # dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified
                colo.res=coloc.abf(D1_l, D2_l)
                colo_res=data.frame(t(colo.res$summary))
                print('|||Coloc result:|||')
                print(colo_res)

                if (colo_res$PP.H4.abf>0.5){

                    colo_df <- data.table(
                      gwas_name = GWAS_name,
                      gwas_input = variant,
                      locus = paste(chromosome1, locus_start, locus_end, sep = '-'),
                      gwas_lead = GWAS_signal,
                      eqtl_name= eQTL_name,
                      gene = qtl1,
                      eqtl_lead = independent_eqtl_SNP_to_contition_on,
                      pp_h4 = colo.res$summary[['PP.H4.abf']],
                      colocolised_snp = subset(colo.res$results, SNP.PP.H4>0.01)$snp
                    )

                    coloc_results <- append(coloc_results, colo_df)
                    jpeg(paste(variant_id, qtl1, chromosome1, '_GWAS_Conditioned_on_', GWAS_signal, '_eQTL_Conditioned_on_', independent_eqtl_SNP_to_contition_on, 'coloc.jpg', sep='_'))
                        sensitivity(colo.res, "H4 > 0.5")
                    dev.off()

                    if(nrow(colo_res)>0){
                    jpeg(paste(variant_id, qtl1, chromosome1, '_GWAS_Conditioned_on_', GWAS_signal, '_eQTL_Conditioned_on_', independent_eqtl_SNP_to_contition_on, 'condiotioned_rplowt.jpg', sep='_'))
                        par(mfrow=c(2,1))
                        coloc::plot_dataset(D1, highlight_list = list(
                          cond = GWAS_signal, coloc = colo_df$colocolised_snp)
                        )
                        par(new=TRUE)
                        print("Outcome")

                        title(variant_id, line = -2, outer = TRUE)
                        coloc::plot_dataset(D2, highlight_list = list(
                          cond = independent_eqtl_SNP_to_contition_on, coloc = colo_df$colocolised_snp)
                        )
                        par(new=TRUE)
                        print("eQTL")
                    dev.off()
                    }
                }
            }
        }
    }
}

coloc_df <- rbindlist(coloc_results)
fwrite(coloc_df, 'coloc_results.csv')
