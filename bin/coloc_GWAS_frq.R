#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GWAS = args[1]
bfile = args[2]
coloc_input_file = args[3]
# GWAS = 'GCST90014122_buildGRCh37.tsv'
gwas_significance_threshold <- 5e-8
eqtl_significance_threshold <- 5e-5

print(GWAS)
print(bfile)

library(data.table)
source('dataIO.R')

return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats = return_list$map
GWAS_name = return_list$GWAS_name
Significant_GWAS_Signals = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$p_value < gwas_significance_threshold]

writeLines(Full_GWAS_Sum_Stats$variant_id, con=paste0(GWAS_name, ".snp.list"))
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
  single_eqtl1 = load_eqtl(val,Full_GWAS_Sum_Stats)
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

plink_cmd <- paste0("plink --bfile ", bfile, "/plink_genotypes --maf 0.0001 --make-bed --freqx --out ", GWAS_name)
rc <- system(plink_cmd)
stopifnot(rc == 0L)