#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
eQTL = args[1] #'samplename'
GWAS =args[2] # '../donor_ids.tsv'
# eQTL="/lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline/work/86/e52af0eb213bf3984783e560bde9c0/Oligodendrocytes.6.gz"
# GWAS="/lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline/work/86/e52af0eb213bf3984783e560bde9c0/GCST90014123_buildGRCh37.tsv"
print('eQTL')
print(eQTL)
print('GWAS')
print(GWAS)
#if(!require("remotes"))
#  install.packages("remotes") # if necessary
#library(remotes)
#install_github("chr1swallace/coloc",build_vignettes=TRUE)
#devtools::install_github("mrcieu/gwasglue")

library(data.table)
library(coloc)
library(susieR)


#### Eqtl data

eqtl.file=eQTL #downloaded from Zenodo
eqtl=fread(eqtl.file)
names(eqtl)=c("gene","SNP","TSS_dist","p","beta")
eqtl$se=abs(eqtl$beta)/sqrt(qchisq(eqtl$p,df = 1,lower.tail = F))
Celltype=tail(strsplit(eqtl.file,split="/")[[1]],n=1)
Celltype=gsub('\\.','_',Celltype)

map=fread(GWAS)
#First gene snp pair
low_eqtl = eqtl[eqtl$p< 5e-3]
all_uq = unique(low_eqtl$gene)
uq_subset = all_uq[10:length(all_uq)]
uq_subset = all_uq
### Get alla available datasets in ieugwasr
# datasets=ieugwasr::gwasinfo()
# datasets=as.data.frame(datasets)
# datasets[grep("Body mass index",datasets$trait),]
print('#####################')
print('#########---We will test:----############')
print(uq_subset)
print('#####################')
print('#####################')

for (qtl1 in uq_subset){
  print(qtl1)
  # qtl1='ZSCAN23_ENSG00000187987'
  single_eqtl=eqtl[eqtl$gene==qtl1,]

  # Add snp info
  
  single_eqtl$ea_allele=map$effect_allele[match(single_eqtl$SNP,map$variant_id)]
  single_eqtl$chromosome=map$chromosome[match(single_eqtl$SNP,map$variant_id)]
  single_eqtl$base_pair_location=map$base_pair_location[match(single_eqtl$SNP,map$variant_id)]
  single_eqtl = single_eqtl[!is.na(single_eqtl$ea_allele)]
  outcome=map[match(single_eqtl$SNP,map$variant_id),]

  if (nrow(single_eqtl)>0){

    ### There is a eQTL and Two possible signals in the BMI GWAS
    #Calulate LD for SNPs in the list

    ld=ieugwasr::ld_matrix(variants = single_eqtl$SNP, bfile = "/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes", plink_bin = "/software/hgi/installs/anaconda3/envs/hgi_base/bin/plink")
    info=as.data.frame(matrix(unlist(apply(t(row.names(ld)),MARGIN = 1,function(x)strsplit(x,split="_"))),ncol=3,byrow = T))
    names(info)=c("SNP","a0","a1")
    #Harmonise datasets to ld matrix
    single_eqtl=single_eqtl[match(info$SNP,single_eqtl$SNP),]
    to.flip=which(single_eqtl$ea_allele==info$a0 & single_eqtl$oth_allele==info$a1)
    if(length(to.flip)>0){
      single_eqtl$beta[to.flip]=single_eqtl$beta[to.flip]*-1
      a1=single_eqtl$ea_allele[to.flip]
      a0=single_eqtl$oth_allele[to.flip]
      single_eqtl$ea_allele[to.flip]=a0
      single_eqtl$oth_allele[to.flip]=a1
    }
    bad=which(single_eqtl$ea_allele!=info$a1 | single_eqtl$oth_allele!=info$a0)
    if(length(bad)>0) single_eqtl=single_eqtl[-bad,] 

    outcome=outcome[match(info$SNP,outcome$variant_id),]
    to.flip=which(outcome$effect_allele==info$a0 & outcome$other_allele==info$a1)
    if( length(to.flip)>0){
      outcome$beta[to.flip]=outcome$beta[to.flip]*-1
      a1=outcome$effect_allele[to.flip]
      a0=outcome$other_allele[to.flip]
      outcome$effect_allele[to.flip]=a0
      outcome$other_allele[to.flip]=a1
    }
    bad=which(outcome$ea!=info$a1 | outcome$nea!=info$a0)
    if(length(bad)>0)  outcome=outcome[-bad,] 

    #Run Susie eqtl
    fitted_eqtl <- susie_rss(z = single_eqtl$beta/single_eqtl$se, R = ld,L = 10,coverage = 0.85)
    print('#####################')
    print('#########---eQTL---############')
    print(summary(fitted_eqtl))
    print('#####################')
    #Run Susie outcome
    fitted_outcome <- susie_rss(z = outcome$beta/outcome$standard_error, R = ld,L = 10,coverage = 0.85)
    print('#####################')
    print('########---Outcome---#############')
    print(summary(fitted_outcome))
    print('#####################')


    if (length(fitted_outcome$sets$cs)>0 && length(fitted_eqtl$sets$cs)>0){
      print('yes')
      # Here create a new folder
      #run susie coloc
      colocalisation=coloc.susie(fitted_eqtl,fitted_outcome)
      if (!is.null(colocalisation$summary)){

        jpeg(paste(qtl1,'_rplot.jpg',sep=''))
        par(mfrow=c(2,1))
        plot(single_eqtl$base_pair_location,-log10(single_eqtl$p),pch=20)
        plot(outcome$base_pair_location,-log10(outcome$p_value),pch=20)
        dev.off()

        # sensitivity(colocalisation,rule = "H4>0.05",plot.manhattans = F)
        write.table(colocalisation, file = paste(qtl1,'_colocalisation.tsv',sep=''),sep='\t')
      }
      
      print('########RESULTS#############')
      print(colocalisation)
      print('#####################')


    }


  }


}
print('Done')
  # outcome=map[list_of_vars]
  # outcome=map[map$variant_id==single_eqtl$SNP,]
  
  
  # rownames(map)<- map$variant_id
  # single_eqtl$oth_allele=map$other_allele[match(single_eqtl$SNP,map$variant_id)]
  # 
  # # single_eqtl$SNP[0:400]
  # outcome=ieugwasr::associations(variants =single_eqtl$SNP[0:500],id =  "ukb-b-1489")
  # outcome=as.data.frame(outcome)
  # 
  # outcome=ieugwasr::associations(variants =single_eqtl$SNP[1501:2000],id =  "ukb-b-1489")
  # outcome2=as.data.frame(outcome)  
  # outcome=ieugwasr::associations(variants =single_eqtl$SNP[2001:3000],id =  "ukb-b-1489")
  # outcome3=as.data.frame(outcome) 
  # outcome = rbind(outcome3,outcome2,outcome1)
#   
#   if(nrow(outcome)>500){
#     outcome=outcome[1:500,]
# }                             #LD is limited to 500 SNPs

# library(qqman)
# jpeg('manh_rplot4.jpg')
# par(mfrow=c(2,1))
# manhattan(outcome, chr='chromosome', bp='base_pair_location', snp='variant_id', p='p_value' )
# manhattan(single_eqtl, chr='chromosome', bp='base_pair_location', snp='SNP', p='p' )
# dev.off()
### make regional p-plot 

# H3 (both traits associated in the locus but independently) highly likely 



