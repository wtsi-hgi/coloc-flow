#!/usr/bin/env Rscript
# rm(list = setdiff(ls(), lsf.str()))
args = commandArgs(trailingOnly=TRUE)
eQTL = args[1] #'samplename'
GWAS =args[2] # '../donor_ids.tsv'
# eQTL="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/eQTLs/all/Inhibitory.neurons.2.gz"
# GWAS="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/Lacular_stroke/GCST90014122_buildGRCh37.tsv"
print('eQTL')
print(eQTL)
print('GWAS')
print(GWAS)
#if(!require("remotes"))
#  install.packages("remotes") # if necessary
#library(remotes)
#install_github("chr1swallace/coloc",build_vignettes=TRUE)
#devtools::install_github("mrcieu/gwasglue")
susie_rss_fit_GWAS <- function(outcome){
    # outcome=variants_of_interest[match(info$SNP,variants_of_interest$variant_id),]
    ld_GWAS=ieugwasr::ld_matrix(variants = outcome$variant_id, bfile = "/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes", plink_bin = "/software/hgi/installs/anaconda3/envs/hgi_base/bin/plink")
    info=as.data.frame(matrix(unlist(apply(t(row.names(ld_GWAS)),MARGIN = 1,function(x)strsplit(x,split="_"))),ncol=3,byrow = T))
    names(info)=c("SNP","a0","a1")
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
    #Run Susie outcome
    fitted_outcome <- susie_rss(z = outcome$beta/outcome$standard_error, R = ld_GWAS,L = 10,coverage = 0.75,max_iter = 150)
    print('#####################')
    print('########---Outcome---#############')
    print(summary(fitted_outcome))
    print('#####################')
  return(fitted_outcome)
}

susie_rss_fit <- function(single_eqtl){

  # eQTL
  ld_eQLT=ieugwasr::ld_matrix(variants = single_eqtl$SNP, bfile = "/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes", plink_bin = "/software/hgi/installs/anaconda3/envs/hgi_base/bin/plink")
  info=as.data.frame(matrix(unlist(apply(t(row.names(ld_eQLT)),MARGIN = 1,function(x)strsplit(x,split="_"))),ncol=3,byrow = T))
  names(info)=c("SNP","a0","a1")
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
  #Run Susie eqtl
  fitted_eqtl <- susie_rss(z = single_eqtl$beta/single_eqtl$se, R = ld_eQLT,L = 10,coverage = 0.75,max_iter = 150)
  # fitted_eqtl <- susie_rss(z = single_eqtl$beta/single_eqtl$se, R = ld_eQLT,L = 50,coverage = 0.75,max_iter = 200)
  print('#####################')
  print('#########---eQTL---############')
  print(summary(fitted_eqtl))
  print('#####################')


  return(fitted_eqtl)
}

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
map2 = map[map$p_value< 5e-8]
#First gene snp pair
# single_eqtl1 = eqtl[eqtl$p< 5e-2] #We do notfilter out the genes.
single_eqtl1 = eqtl


# Here we diver a bit and add the positional info to the SNPs if available.
single_eqtl1$ea_allele=map$effect_allele[match(single_eqtl1$SNP,map$variant_id)]
single_eqtl1$oth_allele=map$other_allele[match(single_eqtl1$SNP,map$variant_id)]
single_eqtl1$chromosome=map$chromosome[match(single_eqtl1$SNP,map$variant_id)]
single_eqtl1$base_pair_location=map$base_pair_location[match(single_eqtl1$SNP,map$variant_id)]
single_eqtl1 = single_eqtl1[!is.na(single_eqtl1$ea_allele)]

# Now we select the variants in a 1 mb window of the gwas hit. 
for (row in 1:nrow(map2)){ 
  base_pair_location <- map2[row, "base_pair_location"]
  variant_id <- map2[row, "variant_id"]$variant_id
  chromosome1 <- map2[row, "chromosome"]$chromosome
  print(row)
  print(paste(variant_id,chromosome1,sep='-'))
  # Select only the eQTLs that are on the same chromosome.
  t3 = single_eqtl1[single_eqtl1$chromosome == chromosome1]
  if(nrow(t3)>0){
      range_min = base_pair_location-1000000
      range_max = base_pair_location+1000000
      variants_of_interest = map[map$chromosome==chromosome1]
      variants_of_interest = variants_of_interest[variants_of_interest$base_pair_location>range_min$base_pair_location & variants_of_interest$base_pair_location<range_max$base_pair_location]
      outcome = variants_of_interest
      fitted_outcome = susie_rss_fit_GWAS(outcome)
      single_eqtl_all=t3[t3$base_pair_location>range_min$base_pair_location & t3$base_pair_location<range_max$base_pair_location]
      all_uq = unique(single_eqtl_all$gene)
      # Check which genes in this region has a signal
      
      # we run the coloc for a gene at a time snce we cant have repeated variants in a ld matrix
      for (qtl1 in all_uq){
        print(qtl1)
        # qtl1='ATP6V1E2_ENSG00000250565'
        # single_eqtl=single_eqtl_all[single_eqtl_all$gene==qtl1,]
        # variants_of_interest3 = single_eqtl[single_eqtl$base_pair_location>range_min$base_pair_location & single_eqtl$base_pair_location<range_max$base_pair_location]
        single_eqtl=t3[t3$gene==qtl1,] 
        # Here loop through each of the individual genes to colocalise their SNPs.
        # all_snps = unique(c(single_eqtl$SNP, variants_of_interest$variant_id))
        # Have to do the mapping for a gene at a time in the locus, since the sushie needs an ld matrix with no repeated SNPs.

        if (min(single_eqtl$p)<5e-4){
                    vars = intersect(outcome$variant_id,single_eqtl$SNP)
          outcome2 = outcome[match(vars,outcome$variant_id),]
          rownames(outcome2) <- outcome2$variant_id
          single_eqtl2 = single_eqtl[match(vars,single_eqtl$SNP),]
          rownames(single_eqtl2) <- single_eqtl2$SNP

          fitted_outcome = susie_rss_fit_GWAS(outcome2)
          fitted_eqtl = susie_rss_fit(single_eqtl2)

          # print(summary(fitted_eqtl))
          # if (length(fitted_eqtl$sets$cs)>0){
              # GWAS
              jpeg(paste('',variant_id,'_',qtl1,'_',chromosome1,'rplowt.jpg',sep=''))
                par(mfrow=c(2,1))
                plot(outcome$base_pair_location,-log10(outcome$p_value),col=ifelse(outcome$base_pair_location %in% c(base_pair_location$base_pair_location), 'red', 'black'),pch=20,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000))
                title(variant_id, line = -2, outer = TRUE)
                plot(single_eqtl$base_pair_location,-log10(single_eqtl$p),col=ifelse(single_eqtl$base_pair_location %in% c(single_eqtl2$base_pair_location), 'black','red'),pch=20,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000))
              dev.off()

              if (length(fitted_outcome$sets$cs)>0 && length(fitted_eqtl$sets$cs)>0){
                print('yes')
                # Here create a new folder
                #run susie coloc
                colocalisation=coloc.susie(fitted_eqtl,fitted_outcome)
                if (!is.null(colocalisation$summary)){
                  
                  write.table(colocalisation, file = paste(qtl1,'_colocalisation.tsv',sep=''),sep='\t')
                  # sensitivity(colocalisation,rule = "H4>0.05",plot.manhattans = F)
                  # write.table(colocalisation, file = paste(qtl1,'_colocalisation2.tsv',sep=''),sep='\t')
                }
                
                print('########RESULTS#############')
                print(colocalisation)
                print('#####################')
              }
          # } #END OF FITTED EQTL
        } #END of CHECKING IF THERE IS A SIGNAL IN EQTL
      } #END OF LOOPING THROUGH QTLS
    } #END OF CHECKIN IF THERE ARE VARIANTS
  } # now select the SNPs that overlap with the 2mb window of GWAS (this may need to be changed to gene selection criteria)


# done
