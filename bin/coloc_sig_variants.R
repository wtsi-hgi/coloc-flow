#!/usr/bin/env Rscript
library(data.table)
library(coloc)
library(susieR)

source('dataIO.R')
source('cojo.R')
source('helpers.R')

# args = commandArgs(trailingOnly=TRUE)
# eQTL = args[1] #'samplename'
# GWAS =args[2] # '../donor_ids.tsv'
# variant = args[3]
eQTL="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/eQTLs/all/Pericytes.17.gz"
eqtl_samples_number = 192
GWAS="GWAS_UKB_logWMHnorm.txt"
variant="rs4588035"
freq_file=paste0(GWAS_name, ".frqx")
bfile='/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes'
# Takes
#1 frq file
#2 variant name
GWAS_name = tools::file_path_sans_ext(basename(GWAS))
eQTL_name = tools::file_path_sans_ext(basename(eQTL))
eQTL_name = gsub("\\.", "_", eQTL_name)

# freqs = read_freqs(freq_file)

return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats = return_list$map
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats)
single_eqtl1 = load_eqtl(eQTL, eqtl_marker_file)

dataset.list=list()
dataset.list$results=list()
coloc_results=c()

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

variant_build <- get_snp_build_version(rs = variant_id, pos = base_pair_location)
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
    D1$N <- unique(Cojo_Dataframe$N)

    for (qtl1 in all_unique_eQTL_signals_in_this_GWAS_range){
        single_eqtl = genes_of_interest[genes_of_interest$gene==qtl1, ]
        single_eqtl2 = single_eqtl
        rownames(single_eqtl2) <- single_eqtl2$SNP
        single_eqtl2$N <- eqtl_samples_number

        if (min(single_eqtl2$p) < eqtl_significance_threshold){
            Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2)

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
                D2$N <- eqtl_samples_number

                SNPs_To_Colocalise = intersect(D1$snp,D2$snp)
                D1_col = D1[D1$snp %in% SNPs_To_Colocalise,]
                D2_col = D2[D2$snp %in% SNPs_To_Colocalise,]

                # dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified
                colo.res=coloc.abf(D1, D2)
                colo_res=data.frame(t(colo.res$summary))
                colo_res$GWAS_hit1=''
                colo_res$eQTL_hit2=''
                colo_res$GWAS_Conditioned_SNP=GWAS_signal
                colo_res$eQTL=qtl1
                colo_res$eQTL_Conditioned_SNP=independent_eqtl_SNP_to_contition_on
                colo_res$GWAS_Range=paste(locus_start, ':', locus_end)
                print('|||Coloc result:|||')
                print(colo_res)

                if (colo_res$PP.H4.abf>0.5){

                    coloc_results=rbind(coloc_results,colo_res)
                    jpeg(paste(variant_id,'_',qtl1,'_',chromosome1,'__GWAS_Conditioned_on__',GWAS_signal,'__eQTL_Conditioned_on__',independent_eqtl_SNP_to_contition_on,'_coloc.jpg',sep=''))
                        sensitivity(colo.res,"H4 > 0.5")
                    dev.off()
                    if(nrow(colo_res)>0){
                    jpeg(paste(variant_id,'_',qtl1,'_',chromosome1,'__GWAS_Conditioned_on__',GWAS_signal,'__eQTL_Conditioned_on__',independent_eqtl_SNP_to_contition_on,'condiotioned_rplowt.jpg',sep=''))
                        par(mfrow=c(2,1))
                        plot(D2_col$position,-log10(D2_col$pvalues),col=ifelse(D2_col$snp %in% c(independent_eqtl_SNP_to_contition_on), 'red', 'black'),pch=10,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(D2_col$pvalues))+1)), lty = 0.3,lwd=1)
                        par(new=TRUE)
                        print("Outcome")
                        title(variant_id, line = -2, outer = TRUE)
                        plot(D1_col$position,-log10(D1_col$pvalues),col=ifelse(D1_col$snp %in% c(GWAS_signal), 'red', 'black'),pch=10,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(D1_col$pvalues))+1)), lty = 0.3,lwd=1)
                        par(new=TRUE)
                        print("eQTL")
                    dev.off()
                    }
                }
            }
        }
    }
}
