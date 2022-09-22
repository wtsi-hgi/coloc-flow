#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GWAS = args[1] #'samplename'
bfile = args[2] #'samplename'
coloc_input_file = args[3]
# GWAS = 'GCST90014122_buildGRCh37.tsv'

# bfile = '.'
print(GWAS)
print(bfile)

library(biomaRt)
library(data.table)
library(coloc)
library(susieR)

load_GWAS <- function(GWAS){

  GWAS_ext = tail(strsplit(GWAS,split="\\.")[[1]],n=1)
  GWAS_name = strsplit(tail(strsplit(GWAS,split="\\/")[[1]],n=1),split="\\.")[[1]][1]

  if(GWAS_ext=='zip'){
    print('zip')
    map=fread(paste('unzip -p ',GWAS,sep=''))

  } else if (GWAS_ext=='gz'){
    map=fread(paste('gunzip -cq ',GWAS,sep=''))
  } else{
    map=fread(GWAS)
  }

  #Gwas col rename
  names(map)[names(map) == 'P'] <- "p_value"
  names(map)[names(map) == 'REF'] <- "effect_allele"
  names(map)[names(map) == 'ALT'] <- "other_allele"
  names(map)[names(map) == 'CHROM'] <- "chromosome"
  names(map)[names(map) == 'ID'] <- "variant_id"
  names(map)[names(map) == 'POS'] <- "base_pair_location"
  names(map)[names(map) == 'BETA'] <- "beta"
  names(map)[names(map) == 'SE'] <- "standard_error"

  names(map)[names(map) == 'P-value'] <- "p_value"
  names(map)[names(map) == 'Allele1'] <- "effect_allele"
  names(map)[names(map) == 'Allele2'] <- "other_allele"
  names(map)[names(map) == 'CHR'] <- "chromosome"
  names(map)[names(map) == 'MarkerName'] <- "variant_id"
  names(map)[names(map) == 'pos'] <- "base_pair_location"
  names(map)[names(map) == 'Effect'] <- "beta"
  names(map)[names(map) == 'StdErr'] <- "standard_error"

  names(map)[names(map) == 'P-value'] <- "p_value"
  names(map)[names(map) == 'Allele1'] <- "effect_allele"
  names(map)[names(map) == 'Allele2'] <- "other_allele"
  names(map)[names(map) == 'CHR'] <- "chromosome"
  names(map)[names(map) == 'MarkerName'] <- "variant_id"
  names(map)[names(map) == 'pos'] <- "base_pair_location"
  names(map)[names(map) == 'Effect'] <- "beta"
  names(map)[names(map) == 'StdErr'] <- "standard_error"
  names(map)[names(map) == 'NStudies'] <- "N"


  return_list <- list("map" = map, "GWAS_name" = GWAS_name)
  return(return_list)

}

load_eqtl <- function(eQTL,Full_GWAS_Sum_Stats){
  #### Eqtl data
  eqtl.file=eQTL #downloaded from Zenodo
  eqtl=fread(eqtl.file)
  names(eqtl)=c("gene","SNP","TSS_dist","p","beta")
  eqtl$se=abs(eqtl$beta)/sqrt(qchisq(eqtl$p,df = 1,lower.tail = F))
  Celltype=tail(strsplit(eqtl.file,split="/")[[1]],n=1)
  Celltype=gsub('\\.','_',Celltype)
  single_eqtl1 = eqtl

  # Here we diver a bit and add the positional info to the SNPs if available.
  single_eqtl1$ea_allele=Full_GWAS_Sum_Stats$effect_allele[match(single_eqtl1$SNP,Full_GWAS_Sum_Stats$variant_id)]
  single_eqtl1$oth_allele=Full_GWAS_Sum_Stats$other_allele[match(single_eqtl1$SNP,Full_GWAS_Sum_Stats$variant_id)]
  single_eqtl1$chromosome=Full_GWAS_Sum_Stats$chromosome[match(single_eqtl1$SNP,Full_GWAS_Sum_Stats$variant_id)]
  single_eqtl1$base_pair_location=Full_GWAS_Sum_Stats$base_pair_location[match(single_eqtl1$SNP,Full_GWAS_Sum_Stats$variant_id)]
  single_eqtl1 = single_eqtl1[!is.na(single_eqtl1$ea_allele)]

  #eQTL col rename
  
  names(single_eqtl1)[names(single_eqtl1) == 'p'] <- "p_value"
  names(single_eqtl1)[names(single_eqtl1) == 'ea_allele'] <- "effect_allele"
  names(single_eqtl1)[names(single_eqtl1) == 'oth_allele'] <- "other_allele"
  names(single_eqtl1)[names(single_eqtl1) == 'SNP'] <- "variant_id"
  names(single_eqtl1)[names(single_eqtl1) == 'base_pair_location'] <- "base_pair_location"
  names(single_eqtl1)[names(single_eqtl1) == 'beta'] <- "beta"
  names(single_eqtl1)[names(single_eqtl1) == 'se'] <- "standard_error"
  return (single_eqtl1)
}

return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats=return_list$map
GWAS_name=return_list$GWAS_name
Significant_GWAS_Signals = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$p_value< 5e-8]

write(Full_GWAS_Sum_Stats$variant_id,ncol=1,file=paste0(GWAS_name,".snp.list"),sep = "\t")
Significant_GWAS_Signals$gwas_name=paste(Significant_GWAS_Signals$variant_id,'--',GWAS_name,sep='')

# Here we should loop through the input file eQTLs for the particular GWAS and replicate the table so many times

write.table(Significant_GWAS_Signals,file=paste0(GWAS_name,".sig_signals.list"),sep = "\t",quote = FALSE,row.names = FALSE)
# GWAS_name='GCST90014122_buildGRCh37'
Significant_GWAS_Signals2 = read.table(paste0(GWAS_name,".sig_signals.list"),sep = "\t",header=TRUE)
coloc_input = read.table(coloc_input_file,header=TRUE)
all_eQTLs_associated_with_this_GWAS = coloc_input[coloc_input$GWAS  %like% GWAS_name,]
Data2 = data.table()

for (val in all_eQTLs_associated_with_this_GWAS$eQTL){
  # Here we redune the computational testing burden of spining up and reading in same file multiple times by prereading the files here and seeing whether there is a signal in the ceirtain file on the particular chromosomes where GWAS signal is present.
  single_eqtl1 = load_eqtl(val,Full_GWAS_Sum_Stats)
  single_eqtl2 = single_eqtl1[single_eqtl1$p_value<5e-5]
  uq1 = unique(single_eqtl2$chromosome)
  un2 = unique(Significant_GWAS_Signals2$chromosome)
  int1 = intersect(un2,uq1)
  if(length(int1)>0){
      # We only bind the eQTL file if both contain a signal on a specific chromosome.
      # furthermore we should only select the variants that are on particular chromosomes for the analysis.
      Significant_GWAS_Signals_new = Significant_GWAS_Signals2[Significant_GWAS_Signals2$chromosome %in% unique(single_eqtl2$chromosome),]
      Significant_GWAS_Signals_new$gwas_name2 = paste(Significant_GWAS_Signals_new$gwas_name,'--',val,sep='')
      Data2 = rbind(Data2, Significant_GWAS_Signals2) 
  }
}

write.table(Data2,file=paste0(GWAS_name,"_all_signals.tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
system(paste("plink --bfile ",bfile,"/plink_genotypes --extract ",GWAS_name,".snp.list --maf 0.0001 --make-bed --freqx --out ",GWAS_name,sep=''))
