#!/usr/bin/env Rscript
# rm(list = setdiff(ls(), lsf.str()))
args = commandArgs(trailingOnly=TRUE)
eQTL = args[1] #'samplename'
GWAS =args[2] # '../donor_ids.tsv'
eQTL="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/eQTLs/all/Oligodendrocytes.3.gz"
GWAS="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/Lacular_stroke/GCST90014123_buildGRCh37.tsv"
bfile = "/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes"
tmp_dir = '/lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/tmp2'

print('eQTL')
print(eQTL)
print('GWAS')
print(GWAS)
library(biomaRt)
library(data.table)
library(coloc)
library(susieR)

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
    D3=''
    D3$LD = ld_GWAS
    D3$beta = outcome$beta
    D3$varbeta = outcome$standard_error
    D3$snp = outcome$variant_id
    D3$position = outcome$base_pair_location
    D3$type = 'cc'

    check_dataset(D3,req="LD")
    fitted_outcome <- susie_rss(z = outcome$beta/outcome$standard_error, R = ld_GWAS,L = 10,coverage = 0.95,max_iter = 50,prior_variance = 150,)
    print('#####################')
    print('########---Outcome---#############')
    print(summary(fitted_outcome))
    print('#####################')
  return(fitted_outcome)
}

make_cojo_df <- function(outcome2){
          outcome2$FREQ=freqs$FreqA1[match(outcome2$variant_id,freqs$SNP)]
          idx=which(outcome2$effect_allele!=freqs$effect_allele)
          outcome2$FREQ[idx]=1-outcome2$FREQ[idx]
          # outcome2$se=sqrt(outcome2$varbeta)
          Cojo_Dataframe=outcome2[,c("variant_id","effect_allele","other_allele","FREQ","beta","standard_error","p_value","N")]
          names(Cojo_Dataframe)=c("SNP" , "A1" ,  "A2"  , "freq", "b"  ,  "se" ,  "p" ,   "N")
          Cojo_Dataframe=Cojo_Dataframe[!is.na(Cojo_Dataframe$freq)]# <- 0
          Cojo_Dataframe=Cojo_Dataframe[!Cojo_Dataframe$N <10]# <- 10 #Cojo doesnt like sample sizes smaller than 10
          Cojo_Dataframe = transform(Cojo_Dataframe, freq = as.numeric(freq), 
               N = as.numeric(N),b = as.numeric(b))
    return (Cojo_Dataframe)
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

  return(map)

}

load_eqtl <- function(eQTL){
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
  rownames(single_eqtl1) <- single_eqtl1$SNP
  names(single_eqtl1)[names(single_eqtl1) == 'p'] <- "p_value"
  names(single_eqtl1)[names(single_eqtl1) == 'ea_allele'] <- "effect_allele"
  names(single_eqtl1)[names(single_eqtl1) == 'oth_allele'] <- "other_allele"
  names(single_eqtl1)[names(single_eqtl1) == 'SNP'] <- "variant_id"
  names(single_eqtl1)[names(single_eqtl1) == 'base_pair_location'] <- "base_pair_location"
  names(single_eqtl1)[names(single_eqtl1) == 'beta'] <- "beta"
  names(single_eqtl1)[names(single_eqtl1) == 'se'] <- "standard_error"
}



Full_GWAS_Sum_Stats = load_GWAS(GWAS)
Significant_GWAS_Signals = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$p_value< 5e-8]

write(Full_GWAS_Sum_Stats$variant_id,ncol=1,file=paste0(tmp_dir,'/',GWAS_name,".snp.list"))
system(paste("plink --bfile ",bfile," --extract ",tmp_dir,'/',GWAS_name,".snp.list --maf 0.0001 --make-bed --freqx --out ",tmp_dir,'/',GWAS_name,sep=''))
freqs=fread(paste0(tmp_dir,'/',GWAS_name,".frqx"))
freqs$FreqA1=(freqs$'C(HOM A1)'*2+freqs$'C(HET)')/(2*(rowSums(freqs[,c("C(HOM A1)", "C(HET)", "C(HOM A2)")])))          
#First gene snp pair
# single_eqtl1 = eqtl[eqtl$p< 5e-2] #We do notfilter out the genes.






#############################
#############################
# Overall pattern of the folowing loops:
# Prior to everything we create a bfile based on a Blueprint, uk10k, 1kg projects thats extracted from gwas signal fed in.
# we load the eQTL file
# we load the GWAs sumary stats file and process columns accodringly
# we calculate the frequencies of SNPs.

# 1) then from the summary stats we take each a GWAS signals with a significance < 5e-8
# 2) then we check if the eqtl file that we have loaded contains eQTLs on the same chromosome.
# 3) we select +/-1mb of the GWAS SNP and select the eQTL signalis in the location
# 4) We perform analysis of the number of independent SNPs in this location for each eQTL and GWAS.
# 5) Then for each independent GWAS signal we condition signal to this SNP
# 6) We do the same for each independent eQTL signal - condition.
# 7) Then we perform coloc.
# 8) If coloc is successful we would record the data.
#############################
#############################


# Now we select the variants in a 1 mb window of the gwas hit. 
for (row in 1:nrow(Significant_GWAS_Signals)){ 
  row=83
  base_pair_location <- Significant_GWAS_Signals[row, "base_pair_location"]
  variant_id <- Significant_GWAS_Signals[row, "variant_id"]$variant_id
  chromosome1 <- Significant_GWAS_Signals[row, "chromosome"]$chromosome
  print(paste(variant_id,chromosome1,sep='-'))
  # Select only the eQTLs that are on the same chromosome.
  t3 = single_eqtl1[single_eqtl1$chromosome == chromosome1]
  if(nrow(t3)>0){
      range_min = base_pair_location-1000000
      range_max = base_pair_location+1000000
      variants_of_interest = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$chromosome==chromosome1]
      variants_of_interest = variants_of_interest[variants_of_interest$base_pair_location>range_min$base_pair_location & variants_of_interest$base_pair_location<range_max$base_pair_location]
      outcome = variants_of_interest

      # fitted_outcome = susie_rss_fit_GWAS(outcome)
      single_eqtl_all=t3[t3$base_pair_location>range_min$base_pair_location & t3$base_pair_location<range_max$base_pair_location]
      all_uq = unique(single_eqtl_all$gene)
      # Check which genes in this region has a signal
      
      # we run the coloc for a gene at a time since we cant have repeated variants in a ld matrix
      for (qtl1 in all_uq){
        print(qtl1)
        qtl1='ULK4_ENSG00000168038'
        # single_eqtl=single_eqtl_all[single_eqtl_all$gene==qtl1,]
        # variants_of_interest3 = single_eqtl[single_eqtl$base_pair_location>range_min$base_pair_location & single_eqtl$base_pair_location<range_max$base_pair_location]
        single_eqtl=t3[t3$gene==qtl1,] 
        # Here loop through each of the individual genes to colocalise their SNPs.
        # all_snps = unique(c(single_eqtl$SNP, variants_of_interest$variant_id))
        # Have to do the mapping for a gene at a time in the locus, since the sushie needs an ld matrix with no repeated SNPs.

        
        # break
        # We intersect the variants so that they are pointinmg to the same variant. 

        vars = intersect(outcome$variant_id,single_eqtl$SNP)
        outcome2 = outcome[match(vars,outcome$variant_id),]
        rownames(outcome2) <- outcome2$variant_id
        single_eqtl2 = single_eqtl[match(vars,single_eqtl$SNP),]

        # names(single_eqtl2)[names(single_eqtl2) == 'NStudies'] <- "N"
        single_eqtl2$N = 200   #Check if there is an actual n number.
        if (min(single_eqtl2$p)<5e-5){
          # break
          # single_eqtl = single_eqtl2
          # outcome = outcome2
          # fitted_outcome = susie_rss_fit_GWAS(outcome2)
          # GWAS_name = 'GCST90014123_buildGRCh37'
          # plink --bfile /lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes --extract /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp3/rs1613759.snp.list --maf 0.0001 --make-bed --freqx --out /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp3/rs1613759
          # system(paste0("plink --bfile ",bfile," --extract ",variant_id,".snp.list --maf 0.0001 --make-bed --freqx --out ",variant_id))
          
          # write(map$variant_id,ncol=1,file=paste0('/lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/tmp/',variant_id,".snp.list"))
          # print(paste("plink --bfile /lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes --extract /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/tmp/",variant_id,".snp.list --maf 0.0001 --make-bed --freqx --out /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/plink/",variant_id,sep=''))
            
          # plink --bfile /lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes --extract /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37.snp.list --maf 0.0001 --make-bed --freqx --out /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37
          
          # Conditioning the GWAS signal for the SNPs of interest.
          Cojo_Dataframe = make_cojo_df(outcome2)
          write(Cojo_Dataframe$SNP,ncol=1,file=paste0(tmp_dir,'/',variant_id,".snp.list"))
          write.table(Cojo_Dataframe,file=paste0(tmp_dir,'/',variant_id,"_sum.txt"),row.names=F,quote=F,sep="\t")
          system(paste("gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/plink2/all --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,".snp.list --cojo-file ",tmp_dir,"/",variant_id,"_sum.txt --cojo-slct --out ",tmp_dir,'/',variant_id,"_step1", sep=""))
          independent_SNPs=fread(paste0(tmp_dir,"/",variant_id,"_step1.jma.cojo"))

          # Conditioning the eQTL signal for the SNPs of interest.
          Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2)
          write(Cojo_Dataframe_eqtl$SNP,ncol=1,file=paste0(tmp_dir,'/',variant_id,"_eqtl.snp.list"))
          write.table(Cojo_Dataframe_eqtl,file=paste0(tmp_dir,'/',variant_id,"_eqtl_sum.txt"),row.names=F,quote=F,sep="\t")
          system(paste("gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/plink2/all --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,"_eqtl.snp.list --cojo-file ",tmp_dir,"/",variant_id,"_eqtl_sum.txt --cojo-slct --out ",tmp_dir,'/',variant_id,"eqtl_step1", sep=""))
          independent_SNPs_eQTL=fread(paste0(tmp_dir,"/",variant_id,"eqtl_step1.jma.cojo"))

         
          # if only 1 then run as per:
          #if more than 1 then condition for each.
          # gcta64 --bfile /Users/mo11/Downloads/tmp2/rs9958650 --cojo-p 1e-4 --extract /Users/mo11/Downloads/tmp2/rs9958650.snp.list --cojo-file /Users/mo11/Downloads/tmp2/rs9958650_sum.txt --cojo-slct --out /Users/mo11/Downloads/tmp2/rs9958650_step1
          
          dataset.list=list()
          dataset.list$independent_SNPs=independent_SNPs
          dataset.list$results=list()

          # now we loop through each of the independent SNP combinations from cojo and eQTL to colocalise the signals.
          coloc_results=c()

          for( i_GWAS in 1:nrow(independent_SNPs)){

              # Write to a temp folder, if the conditioning already exists then load it in without conditioning.
              print(paste('GWAS signal is = ',independent_SNPs[i_GWAS]$SNP))
              print(paste('eQTL signal is = ',i_eQTL))

              write(independent_SNPs$SNP[i_GWAS],ncol=1,file=paste0(tmp_dir,"/",independent_SNPs$SNP[i_GWAS],"_independent.snp"))
              system(paste("gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/plink2/all --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,".snp.list --cojo-file ",tmp_dir,"/",variant_id,"_sum.txt --cojo-cond ",tmp_dir,"/",independent_SNPs$SNP[i_GWAS],"_independent.snp --out ",tmp_dir,'/',variant_id,"_step2", sep=""))
              conditioned_dataset=fread(paste0(tmp_dir,'/',variant_id,"_step2.cma.cojo"))
              D1=conditioned_dataset[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
              names(D1)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
              D1$type="quant"
              D1$varbeta=D1$varbeta^2
              D1=na.omit(D1)


            for( i_eQTL in 1:nrow(independent_SNPs_eQTL)){


              write(independent_SNPs_eQTL$SNP[i_eQTL],ncol=1,file=paste0(tmp_dir,"/",independent_SNPs_eQTL$SNP[i_eQTL],"_eqtl_independent.snp"))
              system(paste("gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/plink2/all --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,"_eqtl.snp.list --cojo-file ",tmp_dir,"/",variant_id,"_eqtl_sum.txt --cojo-cond ",tmp_dir,"/",independent_SNPs_eQTL$SNP[i_eQTL],"_eqtl_independent.snp --out ",tmp_dir,'/',variant_id,"eqtl_step2", sep=""))
              conditioned_dataset_eQTL=fread(paste0(tmp_dir,'/',variant_id,"eqtl_step2.cma.cojo"))
              D2=conditioned_dataset_eQTL[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
              names(D2)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
              D2$type="quant"
              D2$varbeta=D2$varbeta^2
              D2=na.omit(D2)
              colo.res=coloc.abf(D1,D2)
              

              jpeg(paste(tmp_dir,'/',variant_id,'_',qtl1,'_',chromosome1,'condiotioned_rplowt.jpg',sep=''))
                par(mfrow=c(2,1))
                plot(conditioned_dataset$bp,-log10(conditioned_dataset$p),col=ifelse(conditioned_dataset$bp %in% c(base_pair_location$base_pair_location), 'red', 'black'),pch=19,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(conditioned_dataset$p))+1)), lty = 1,lwd=1)
                par(new=TRUE)
                print("Outcome")
                title(variant_id, line = -2, outer = TRUE)
                plot(conditioned_dataset_eQTL$bp,-log10(conditioned_dataset_eQTL$p),col=ifelse(conditioned_dataset_eQTL$bp %in% c(base_pair_location$base_pair_location), 'red', 'black'),pch=19,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(conditioned_dataset_eQTL$p))+1)), lty = 1,lwd=1)
                par(new=TRUE)
                print("eQTL")
              dev.off()

              colo.res=data.frame(t(colo.res$summary))
              colo.res$hit1=independent_SNPs[i_GWAS]$SNP
              colo.res$hit2=independent_SNPs_eQTL[i_eQTL]$SNP
              coloc_results=rbind(coloc_results,colo.res)
            }
          }


        }
          # gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37  --extract /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37.snp.list  --cojo-file /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37_sum.txt  --cojo-cond 455666_independent.snp --out 455666_step2
          # gcta --bfile /lustre/scratch123/hgi/projects/bhf_finemap/coloc/tmp2/GCST90014123_buildGRCh37 --cojo-p 1e-4 --extract 455666.snp.list  --cojo-file 455666_sum.txt --cojo-slct --cojo-cond 455666_independent.snp --out 455666_step2
                      

          fitted_eqtl = susie_rss_fit(single_eqtl2)

          # print(summary(fitted_eqtl))
          # print(summary(fitted_outcome))
          # if (length(fitted_eqtl$sets$cs)>0){
              # GWAS
              jpeg(paste(tmp_dir,'/',variant_id,'_',qtl1,'_',chromosome1,'rplowt.jpg',sep=''))
                par(mfrow=c(2,1))
                plot(outcome2$base_pair_location,-log10(outcome2$p_value),col=ifelse(outcome2$base_pair_location %in% c(base_pair_location$base_pair_location), 'red', 'black'),pch=19,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(outcome$p_value))+1)), lty = 1,lwd=1)
                par(new=TRUE)
                # plot(outcome[2319]$base_pair_location,-log10(outcome[2319]$p_value),col='green',, cex=1.5,lty = 2,pch=5,lwd=1,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(outcome$p_value))+1)),ylab="",yaxt="n",xlab="",xaxt="n")
                print("Outcome")
                # for (it in 1:length(fitted_outcome$sets$cs)){
                #   r = ls(fitted_outcome$sets$cs[it])[1]
                #   # print(r)
                #   # fitted_eqtl$sets$cs[it][[r]]  
                #   for (t in fitted_outcome$sets$cs[it][[r]]){
                #     # print(single_eqtl[t]$base_pair_location)
                #     print(outcome2[t]$variant_id)
                #     plot(outcome2[t]$base_pair_location,-log10(outcome2[t]$p),col='green',cex=1.5,lty = 2,pch=5,lwd=1,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(outcome2$p_value))+1)),ylab="",yaxt="n",xlab="",xaxt="n")
                #     par(new=TRUE)
                #   }
                # }
                title(variant_id, line = -2, outer = TRUE)
                
                plot(single_eqtl2$base_pair_location,-log10(single_eqtl2$p),col=ifelse(single_eqtl2$base_pair_location %in% c(single_eqtl2$base_pair_location), 'black','red'),pch=20,ylim=c(0,round(max(-log10(single_eqtl2$p))+1)),xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000))
                par(new=TRUE)
                print("eQTL")
                # for (it in 1:length(fitted_eqtl$sets$cs)){
                #   r = ls(fitted_eqtl$sets$cs[it])[1]
                #   # print(r)
                #   # fitted_eqtl$sets$cs[it][[r]]  
                #   for (t in fitted_eqtl$sets$cs[it][[r]]){
                #     # print(single_eqtl[t]$base_pair_location)
                #     print(single_eqtl2[t]$SNP)
                #     plot(single_eqtl2[t]$base_pair_location,-log10(single_eqtl2[t]$p),col='green',cex=1.5,lty = 2,pch=5,lwd=1,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(single_eqtl2$p))+1)),ylab="",yaxt="n",xlab="",xaxt="n")
                #     par(new=TRUE)
                #   }
                # }
              dev.off()

              jpeg(paste(tmp_dir,'/',variant_id,'_',qtl1,'_',chromosome1,'cond_rplowt.jpg',sep=''))
                par(mfrow=c(2,1))
                plot(step2.res$bp,-log10(step2.res$p),col=ifelse(step2.res$bp %in% c(base_pair_location$base_pair_location), 'red', 'black'),pch=19,xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000),ylim=c(0,round(max(-log10(outcome$p_value))+1)), lty = 1,lwd=1)
                par(new=TRUE)
                print("Outcome")
                title(variant_id, line = -2, outer = TRUE)
                plot(single_eqtl2$base_pair_location,-log10(single_eqtl2$p),col=ifelse(single_eqtl2$base_pair_location %in% c(single_eqtl2$base_pair_location), 'black','red'),pch=20,ylim=c(0,round(max(-log10(single_eqtl2$p))+1)),xlim = c(range_min$base_pair_location-1000000,range_max$base_pair_location+1000000))
                par(new=TRUE)
                print("eQTL")

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
      } #END OF LOOPING THROUGH QTLS - 
    } #END OF CHECKIN IF THERE ARE VARIANTS
  } # now select the SNPs that overlap with the 2mb window of GWAS (this may need to be changed to gene selection criteria)
# done
