#!/usr/bin/env Rscript
library(optparse)
library(data.table)

option_list <- list(
    make_option('--gwas', action="store", help="path to GWAS summary statistic"),
    make_option('--input', action="store", help="path to input file"),
    make_option('--eqtl_snps', action="store", help = "path to eqtl snp_pos.txt file")
)
args <- parse_args(OptionParser(option_list=option_list))

GWAS = args$gwas
coloc_input_file = args$input
# GWAS = 'GCST90014122_buildGRCh37.tsv'

print(GWAS)

source('dataIO.R')
source('helpers.R')

return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats = return_list$map
GWAS_name = return_list$GWAS_name
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats)

Significant_GWAS_Signals$gwas_name = paste(Significant_GWAS_Signals$variant_id, GWAS, sep='--')

# Here we should loop through the input file eQTLs for the particular GWAS and replicate the table so many times

fwrite(Significant_GWAS_Signals, file=paste0(GWAS_name,".sig_signals.list"), sep = "\t", quote = FALSE, row.names = FALSE)
# GWAS_name='GCST90014122_buildGRCh37'
Significant_GWAS_Signals2 = copy(Significant_GWAS_Signals)
coloc_input = fread(coloc_input_file, header=TRUE)
all_eQTLs_associated_with_this_GWAS = coloc_input[coloc_input$GWAS %like% GWAS_name, ]

data_list <- lapply(all_eQTLs_associated_with_this_GWAS$eQTL, function(val){
  # Here we reduce the computational testing burden of spining up and reading in same file multiple times
  # by prereading the files here and seeing whether there is a signal in the ceirtain file on the particular chromosomes where GWAS signal is present.
  single_eqtl1 = load_eqtl(val, args$eqtl_snps)
  single_eqtl2 = single_eqtl1[single_eqtl1$p_value < eqtl_significance_threshold]
  uq1 = unique(single_eqtl2$chromosome)
  un2 = unique(Significant_GWAS_Signals2$chromosome)
  int1 = intersect(un2, uq1)
  if(length(int1)>0){
      # We only bind the eQTL file if both contain a signal on a specific chromosome.
      # furthermore we should only select the variants that are on particular chromosomes for the analysis.
      Significant_GWAS_Signals_new = Significant_GWAS_Signals2[Significant_GWAS_Signals2$chromosome %in% uq1, ]
      Significant_GWAS_Signals_new$gwas_name2 = paste(Significant_GWAS_Signals_new$gwas_name, val, sep='--')
  } else {
      Significant_GWAS_Signals_new <- data.table()
  }
  return(Significant_GWAS_Signals_new)
})

Data2 <- rbindlist(data_list)

fwrite(Data2, file=paste0(GWAS_name,"_all_signals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
