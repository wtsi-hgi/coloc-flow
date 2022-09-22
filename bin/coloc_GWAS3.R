#!/usr/bin/env Rscript
# rm(list = setdiff(ls(), lsf.str()))
# args = commandArgs(trailingOnly=TRUE)
# eQTL = args[1] #'samplename'
# GWAS =args[2] # '../donor_ids.tsv'
eQTL="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/eQTLs/all/Oligodendrocytes.3.gz"
GWAS="/lustre/scratch123/hgi/projects/bhf_finemap/summary_stats/Lacular_stroke/GCST90014123_buildGRCh37.tsv"
bfile = "/lustre/scratch123/hgi/projects/bhf_finemap/imputation/uk10k_1000g_blueprint/plink_genotypes/plink_genotypes"
tmp_dir = '/lustre/scratch123/hgi/projects/bhf_finemap/coloc/bin/tmp6'

print(eQTL)
print(GWAS)

library(biomaRt)
library(data.table)
library(coloc)
library(susieR)

susie_rss_fit_GWAS <- function(outcome,bfile){
    # outcome=variants_of_interest[match(info$SNP,variants_of_interest$variant_id),]
    ld_GWAS=ieugwasr::ld_matrix(variants = outcome$variant_id, bfile = bfile, plink_bin = "/software/hgi/installs/anaconda3/envs/hgi_base/bin/plink")
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

make_cojo_df <- function(GWAS_matched_SNPS_with_eQTL){
          GWAS_matched_SNPS_with_eQTL$FREQ=freqs$FreqA1[match(GWAS_matched_SNPS_with_eQTL$variant_id,freqs$SNP)]
          idx=which(GWAS_matched_SNPS_with_eQTL$effect_allele!=freqs$effect_allele)
          GWAS_matched_SNPS_with_eQTL$FREQ[idx]=1-GWAS_matched_SNPS_with_eQTL$FREQ[idx]
          # GWAS_matched_SNPS_with_eQTL$se=sqrt(GWAS_matched_SNPS_with_eQTL$varbeta)
          Cojo_Dataframe=GWAS_matched_SNPS_with_eQTL[,c("variant_id","effect_allele","other_allele","FREQ","beta","standard_error","p_value","N")]
          names(Cojo_Dataframe)=c("SNP" , "A1" ,  "A2"  , "freq", "b"  ,  "se" ,  "p" ,   "N")
          Cojo_Dataframe=Cojo_Dataframe[!is.na(Cojo_Dataframe$freq)]# <- 0
          Cojo_Dataframe=Cojo_Dataframe[!Cojo_Dataframe$N <10]# <- 10 #Cojo doesnt like sample sizes smaller than 10
          Cojo_Dataframe = transform(Cojo_Dataframe, freq = as.numeric(freq), N = as.numeric(N),b = as.numeric(b))
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

eQTL_name = strsplit(tail(strsplit(eQTL,split="\\/")[[1]],n=1),split="\\.")[[1]]
eQTL_name = paste(eQTL_name[-length(eQTL_name)],collapse = '_')
return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats=return_list$map
GWAS_name=return_list$GWAS_name
print('Here!!3')
single_eqtl1 = load_eqtl(eQTL,Full_GWAS_Sum_Stats)
print('Here!!32')
Significant_GWAS_Signals = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$p_value< 5e-8]
print('Here!!33')
write(Full_GWAS_Sum_Stats$variant_id,ncol=1,file=paste0(tmp_dir,'/',GWAS_name,".snp.list"))
system(paste("plink --bfile ",bfile," --extract ",tmp_dir,'/',GWAS_name,".snp.list --maf 0.0001 --make-bed --freqx --out ",tmp_dir,'/',GWAS_name,sep=''))

# Frequencies for the same GWAS sum stats may be used in multiple places, hence have to make this a sperate nf process.
freqs=fread(paste0(tmp_dir,'/',GWAS_name,".frqx"))
freqs$FreqA1=(freqs$'C(HOM A1)'*2+freqs$'C(HET)')/(2*(rowSums(freqs[,c("C(HOM A1)", "C(HET)", "C(HOM A2)")])))          

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
print('Here!!2')
dataset.list=list()
dataset.list$results=list()
coloc_results=c()
print('Here!!22')
# Now we select the variants in a 1 mb window of the gwas hit. 
for (row in 1:nrow(Significant_GWAS_Signals)){ 
  # row=83
  print(row)
  print('Here!!')
  base_pair_location <- Significant_GWAS_Signals[row, "base_pair_location"]
  variant_id <- Significant_GWAS_Signals[row, "variant_id"]$variant_id
  chromosome1 <- Significant_GWAS_Signals[row, "chromosome"]$chromosome
  print(paste('Running GWAS variant ',variant_id))
  # Select only the eQTLs that are on the same chromosome.
  eQTL_singals_on_GWAS_SNP_chtomosome = single_eqtl1[single_eqtl1$chromosome == chromosome1]
  if(nrow(eQTL_singals_on_GWAS_SNP_chtomosome)>0){
        range_min = base_pair_location-1000000
        range_max = base_pair_location+1000000
        variants_of_interest = Full_GWAS_Sum_Stats[Full_GWAS_Sum_Stats$chromosome==chromosome1]
        variants_of_interest = variants_of_interest[variants_of_interest$base_pair_location>range_min$base_pair_location & variants_of_interest$base_pair_location<range_max$base_pair_location]
        outcome = variants_of_interest
        # Check for independent signals in this range GWAS
        # question here - do we want to condition a full GWAS range, or only the SNPs that overlap?
        Cojo_Dataframe = make_cojo_df(outcome)
        write(Cojo_Dataframe$SNP,ncol=1,file=paste0(tmp_dir,'/',variant_id,'_',GWAS_name,".snp.list"))
        write.table(Cojo_Dataframe,file=paste0(tmp_dir,'/',variant_id,'_',GWAS_name,"_sum.txt"),row.names=F,quote=F,sep="\t")
        system(paste("gcta --bfile ",tmp_dir,'/',GWAS_name," --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,'_',GWAS_name,".snp.list --cojo-file ",tmp_dir,"/",variant_id,'_',GWAS_name,"_sum.txt --cojo-slct --out ",tmp_dir,'/',variant_id,'_',GWAS_name,"_step1", sep=""))
        independent_SNPs=fread(paste0(tmp_dir,"/",variant_id,'_',GWAS_name,"_step1.jma.cojo"))

        for( i_GWAS in 1:nrow(independent_SNPs)){

            print(paste('GWAS signal is = ',variant_id))
            # i_GWAS=1
            GWAS_signal = independent_SNPs[i_GWAS]$SNP
            print(i_GWAS)
            print(paste('Conditioning the gwas range to ',GWAS_signal))
        
            write(GWAS_signal,ncol=1,file=paste0(tmp_dir,"/",GWAS_signal,"_independent.snp"))
            system(paste("gcta --bfile ",tmp_dir,'/',GWAS_name," --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,'_',GWAS_name,".snp.list --cojo-file ",tmp_dir,"/",variant_id,'_',GWAS_name,"_sum.txt --cojo-cond ",tmp_dir,"/",GWAS_signal,"_independent.snp --out ",tmp_dir,'/',variant_id,'_',GWAS_name,"_step2", sep=""))
            variant_id = 'rs9837273'
            GWAS_name = 'GCST90014123_buildGRCh37'
            GWAS_signal ='rs9831404'
            conditioned_dataset=fread(paste0(tmp_dir,'/',variant_id,'_',GWAS_name,"_step2.cma.cojo"))
            # var id = rs9837273
            # GCST90014123_buildGRCh37

            # Conditioned dataset doesnt nontain the SNP that we condition the data to, this needs to be included. 
            conditioned_dataset_condSNP=fread(paste0(tmp_dir,'/',variant_id,'_',GWAS_name,"_step1.jma.cojo"))
            conditioned_dataset_condSNP=conditioned_dataset_condSNP[conditioned_dataset_condSNP$SNP == GWAS_signal]
                        
            conditioned_dataset_condSNP = conditioned_dataset_condSNP[, !c("LD_r","pJ","bJ","bJ_se")]
            conditioned_dataset_condSNP$pC = conditioned_dataset_condSNP$p
            conditioned_dataset_condSNP$bC = conditioned_dataset_condSNP$b
            conditioned_dataset_condSNP$bC_se = conditioned_dataset_condSNP$se
            conditioned_dataset = rbind(conditioned_dataset, conditioned_dataset_condSNP) 
            # some SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample, hence are removed.
            D1=conditioned_dataset[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
            names(D1)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
            D1$type="quant"
            D1$varbeta=D1$varbeta^2
            D1=na.omit(D1)
        
            single_eqtl_all=eQTL_singals_on_GWAS_SNP_chtomosome[eQTL_singals_on_GWAS_SNP_chtomosome$base_pair_location>range_min$base_pair_location & eQTL_singals_on_GWAS_SNP_chtomosome$base_pair_location<range_max$base_pair_location]
            all_unique_eQTL_signals_in_this_GWAS_range = unique(single_eqtl_all$gene)

            for (qtl1 in all_unique_eQTL_signals_in_this_GWAS_range){
                print(qtl1)
                # qtl1='ABHD5_ENSG00000011198'
                single_eqtl=eQTL_singals_on_GWAS_SNP_chtomosome[eQTL_singals_on_GWAS_SNP_chtomosome$gene==qtl1,] 
                single_eqtl2 = single_eqtl
                rownames(single_eqtl2) <- single_eqtl2$SNP
                single_eqtl2$N = 200   #Check if there is an actual n number.
                if (min(single_eqtl2$p)<5e-5){
                    Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2)
                    write(Cojo_Dataframe_eqtl$SNP,ncol=1,file=paste0(tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list"))
                    write.table(Cojo_Dataframe_eqtl,file=paste0(tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt"),row.names=F,quote=F,sep="\t")
                    system(paste("gcta --bfile ",tmp_dir,'/',GWAS_name," --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list --cojo-file ",tmp_dir,"/",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt --cojo-slct --out ",tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1", sep=""))
                    independent_SNPs_eQTL=fread(paste0(tmp_dir,"/",variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1.jma.cojo"))
                    for( i_eQTL in 1:nrow(independent_SNPs_eQTL)){
                        print(i_eQTL)
                        # i_eQTL=1
                        # issue that there is no snp that we conditioned this on.
                        independent_eqtl_SNP_to_contition_on = independent_SNPs_eQTL[i_eQTL]$SNP
                        write(independent_eqtl_SNP_to_contition_on,ncol=1,file=paste0(tmp_dir,"/",independent_eqtl_SNP_to_contition_on,"_eqtl_independent.snp"))
                        system(paste("gcta --bfile ",tmp_dir,'/',GWAS_name," --cojo-p 1e-4 --extract ",tmp_dir,"/",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list --cojo-file ",tmp_dir,"/",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt --cojo-cond ",tmp_dir,"/",independent_eqtl_SNP_to_contition_on,"_eqtl_independent.snp --out ",tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"_",independent_eqtl_SNP_to_contition_on,"eqtl_step2", sep=""))
                        
                        qtl1='ULK4_ENSG00000168038'
                        eQTL_name='Oligodendrocytes_3'
                        independent_eqtl_SNP_to_contition_on='rs2371592'
                        chromosome1='3'
                        conditioned_dataset_eQTL=fread(paste0(tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"_",independent_eqtl_SNP_to_contition_on,"eqtl_step2.cma.cojo"))
                        
                        # Conditioned dataset doesnt nontain the SNP that we condition the data to, this needs to be included. 
                        conditioned_dataset_eQTL_condSNP=fread(paste0(tmp_dir,'/',variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1.jma.cojo"))
                        conditioned_dataset_eQTL_condSNP=conditioned_dataset_eQTL_condSNP[conditioned_dataset_eQTL_condSNP$SNP == independent_eqtl_SNP_to_contition_on]
                        conditioned_dataset_eQTL_condSNP = conditioned_dataset_eQTL_condSNP[, !c("LD_r","pJ","bJ","bJ_se")]
                        conditioned_dataset_eQTL_condSNP$pC = conditioned_dataset_eQTL_condSNP$p
                        conditioned_dataset_eQTL_condSNP$bC = conditioned_dataset_eQTL_condSNP$b
                        conditioned_dataset_eQTL_condSNP$bC_se = conditioned_dataset_eQTL_condSNP$se
                        conditioned_dataset_eQTL = rbind(conditioned_dataset_eQTL, conditioned_dataset_eQTL_condSNP) 
                        
                        
                        D2=conditioned_dataset_eQTL[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
                        names(D2)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
                        D2$type="quant"
                        D2$varbeta=D2$varbeta^2
                        D2=na.omit(D2)
                        SNPs_To_Colocalise = intersect(D1$snp,D2$snp)
                        # D1_col = D1$snp == SNPs_To_Colocalise
                        D1_col =D1[D1$snp %in% SNPs_To_Colocalise,]
                        D2_col =D2[D2$snp %in% SNPs_To_Colocalise,]

                        # dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified
                        colo.res=coloc.abf(D1,D2)
                        colo_res=data.frame(t(colo.res$summary))
                        colo_res$GWAS_hit1=''
                        colo_res$eQTL_hit2=''
                        colo_res$GWAS_Conditioned_SNP=GWAS_signal
                        colo_res$eQTL=qtl1
                        colo_res$eQTL_Conditioned_SNP=independent_eqtl_SNP_to_contition_on
                        colo_res$GWAS_Range=paste(range_min,':',range_max)
                        print('|||Coloc result:|||')
                        print(colo_res)
                        if (colo_res$PP.H4.abf>0.5){

                        
                          coloc_results=rbind(coloc_results,colo_res)
                          jpeg(paste(tmp_dir,'/',variant_id,'_',qtl1,'_',chromosome1,'__GWAS_Conditioned_on__',GWAS_signal,'__eQTL_Conditioned_on__',independent_eqtl_SNP_to_contition_on,'_coloc.jpg',sep=''))
                                sensitivity(colo.res,"H4 > 0.5")
                          dev.off()

                          if(nrow(colo_res)>0){

                            jpeg(paste(tmp_dir,'/',variant_id,'_',qtl1,'_',chromosome1,'__GWAS_Conditioned_on__',GWAS_signal,'__eQTL_Conditioned_on__',independent_eqtl_SNP_to_contition_on,'condiotioned_rplowt.jpg',sep=''))
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
    }
}







