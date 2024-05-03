#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(dplyr)
library(data.table)
requireNamespace('dplyr')
requireNamespace('tidyr')
library("stringr") 
eqtl_significance_threshold <- 5e-5
gwas_significance_threshold <- 5e-7
option_list <- list(
    make_option('--gwas', action="store", help="path to GWAS summary statistic"),
    make_option('--eqtl_fofn', action="store", help="path to fofn of eqtls"),
    make_option('--eqtl_snps', action="store", help = "path to eqtl snp_pos.txt file")
)
args <- parse_args(OptionParser(option_list=option_list))
# args$gwas = 'GIGASTROKE_SVS_EUR_hg19_harmonised.tsv.gz'
# args$eqtl_fofn = 'eqtl.12.list'
# args$eqtl_snps = 'snp_pos.txt'
GWAS = fs::link_path(args$gwas)
# print(GWAS)
# GWAS = fs::link_path('GIGASTROKE_LAS_EUR_hg19_harmonised.tsv.gz')
source('dataIO.R')
source('helpers.R')
eqtl_marker_data <- read_eqtl_marker_file(args$eqtl_snps, build = 'hg38') # Dependant on which build we are using we preload the SNP variant file to ease the mapping between different builds.

Full_GWAS_Sum_Stats = load_eqtl(GWAS, marker.data = eqtl_marker_data, build = 'hg19')
GWAS_name = tools::file_path_sans_ext(basename(GWAS), compression = T)
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats,threshold=gwas_significance_threshold)
groups <- make_gwas_groups(Significant_GWAS_Signals, window = 1e6) # groups the significant GWAS signals in groups accoring to vindo to reduce the testing burden.

# Here we should loop through the input file eQTLs for the particular GWAS and replicate the table so many times

fwrite(Significant_GWAS_Signals, file=paste0(GWAS_name,".sig_signals.list"), sep = "\t", quote = FALSE, row.names = FALSE)
Significant_GWAS_Signals2 = copy(Significant_GWAS_Signals)

eqtls <- readLines(args$eqtl_fofn)
# eqtls <- readLines('eqtl.14.list' )


data_list <- lapply(eqtls, function(val){
  print(val)
  eqtl.file = val
  # val='/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/eQTLs/all/Oligodendrocytes.12.gz'
  # val='/scratch/cellfunc/shared/HUVEC_RNAseq/eQTLs_norm_counts/TensorQTL_eQTLS/general/nom_output/cis_nominal1.cis_qtl_pairs.chr1.tsv'
  # Here we reduce the computational testing burden of spining up and reading in same file multiple times
  # by prereading the files here and seeing whether there is a signal in the ceirtain file on the particular chromosomes where GWAS signal is present.
  single_eqtl2 = load_eqtl(eqtl.file, marker.data = eqtl_marker_data, build = 'hg38')
  single_eqtl2 <- get_gwas_significant_signals(single_eqtl2,threshold=eqtl_significance_threshold)
  
  uq1 = unique(single_eqtl2$chromosome)
  un2 = unique(Significant_GWAS_Signals2$chromosome)
  int1 = intersect(un2, uq1)
  if(length(int1)>0){
      # break
      # We only bind the eQTL file if both contain a signal on a specific chromosome.
      # furthermore we should only select the variants that are on particular chromosomes for the analysis.
      Significant_GWAS_Signals_new = Significant_GWAS_Signals2[Significant_GWAS_Signals2$chromosome %in% uq1, ]
      Significant_GWAS_Signals_new$base_pair_location_end = Significant_GWAS_Signals_new$base_pair_location

      # from each window of significant markers choose the middle one
      data.table::foverlaps(Significant_GWAS_Signals_new, groups,
                      type = 'within',
                      by.x = c('chromosome', 'base_pair_location', 'base_pair_location_end')) %>%
        add_count(SNP) %>%
        filter(n == 1) -> group_markers

      if(nrow(group_markers) == 0){
        stop('No biallelic markers in significance region!')
      }

      group_markers %>%
          group_by(group_id) %>%
          arrange(base_pair_location) %>%
          summarise(SNP = SNP[ceiling(n()/2)]) -> representative_snps
      representative_snps$gwas_name = paste(representative_snps$SNP, GWAS, sep='--')
      representative_snps$gwas_name2 = paste(representative_snps$gwas_name, val, sep='--')
  } else {
      representative_snps <- data.table()
  }
  return(representative_snps)
})

Data2 <- rbindlist(data_list)

if(ncol(Data2) == 0) Data2 <- data.table(gwas_name2 = character())

fwrite(Data2, file=paste0(GWAS_name, "_all_signals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
