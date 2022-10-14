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
GWAS="GWAS_UKB_logWMHnorm.txt"
variant="rs4588035"
freq_file=paste0(GWAS_name, ".frqx")
# Takes
#1 frq file
#2 variant name
GWAS_name = tools::file_path_sans_ext(basename(GWAS))
eQTL_name = tools::file_path_sans_ext(basename(eQTL))
eQTL_name = gsub("\\.", "_", eQTL_name)

freqs = read_freqs(freq_file)

return_list = load_GWAS(GWAS)
Full_GWAS_Sum_Stats = return_list$map
Significant_GWAS_Signals <- get_gwas_significant_signals(Full_GWAS_Sum_Stats)
single_eqtl1 = load_eqtl(eQTL, Full_GWAS_Sum_Stats)

dataset.list=list()
dataset.list$results=list()
coloc_results=c()

row1 = Significant_GWAS_Signals[Significant_GWAS_Signals$variant_id == variant,]
base_pair_location <- row1[["base_pair_location"]]
variant_id <- row1[["variant_id"]]
chromosome1 <- row1[["chromosome"]]
print(paste('Running GWAS variant', variant_id))

range_min <- base_pair_location - 1000000
range_max <- base_pair_location + 1000000
variants_of_interest <- dplyr::filter(Full_GWAS_Sum_Stats,
    chromosome == chromosome1,
    base_pair_location > range_min,
    base_pair_location < range_max
)

# GWAS_matched_SNPS_with_eQTL=variants_of_interest
Cojo_Dataframe <- make_cojo_df(variants_of_interest, freqs = freqs)

marker_filename <- paste0(variant_id, '_', GWAS_name, ".snp.list")
writeLines(Cojo_Dataframe$SNP, con=marker_filename)

cojo_filename <- paste0(variant_id, '_', GWAS_name, "_sum.txt")
fwrite(Cojo_Dataframe, file=cojo_filename, row.names=F, quote=F, sep="\t")

cojo_out <- run_cojo(
    bfile = paste0("Filtered_", GWAS_name, "/Filtered_", GWAS_name),
    marker_list = marker_filename,
    summary_stat = cojo_filename,
    out_prefix = paste0(variant_id, '_', GWAS_name,"_step1")
)

independent_SNPs = fread(cojo_out)

for( i_GWAS in 1:nrow(independent_SNPs)){
    GWAS_signal = independent_SNPs[i_GWAS]$SNP
    write(GWAS_signal,ncol=1,file=paste0(GWAS_signal,"_independent.snp"))
    system(paste("gcta --bfile ",GWAS_name," --cojo-p 1e-4 --extract ",variant_id,'_',GWAS_name,".snp.list --cojo-file ",variant_id,'_',GWAS_name,"_sum.txt --cojo-cond ",GWAS_signal,"_independent.snp --out ",variant_id,'_',GWAS_name,"_step2", sep=""))
    conditioned_dataset=fread(paste0(variant_id,'_',GWAS_name,"_step2.cma.cojo"))
    conditioned_dataset_condSNP=fread(paste0(variant_id,'_',GWAS_name,"_step1.jma.cojo"))
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
        single_eqtl=eQTL_singals_on_GWAS_SNP_chtomosome[eQTL_singals_on_GWAS_SNP_chtomosome$gene==qtl1,]
        single_eqtl2 = single_eqtl
        rownames(single_eqtl2) <- single_eqtl2$SNP
        single_eqtl2$N = 200   #Check if there is an actual n number.

        if (min(single_eqtl2$p)<5e-5){
            Cojo_Dataframe_eqtl = make_cojo_df(single_eqtl2)
            write(Cojo_Dataframe_eqtl$SNP,ncol=1,file=paste0(variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list"))
            write.table(Cojo_Dataframe_eqtl,file=paste0(variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt"),row.names=F,quote=F,sep="\t")
            system(paste("gcta --bfile ",GWAS_name," --cojo-p 1e-4 --extract ",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list --cojo-file ",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt --cojo-slct --out ",variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1", sep=""))
            independent_SNPs_eQTL=fread(paste0(variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1.jma.cojo"))
            for( i_eQTL in 1:nrow(independent_SNPs_eQTL)){
                print(i_eQTL)
                independent_eqtl_SNP_to_contition_on = independent_SNPs_eQTL[i_eQTL]$SNP
                write(independent_eqtl_SNP_to_contition_on,ncol=1,file=paste0(independent_eqtl_SNP_to_contition_on,"_eqtl_independent.snp"))
                system(paste("gcta --bfile ",GWAS_name," --cojo-p 1e-4 --extract ",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl.snp.list --cojo-file ",variant_id,"_",qtl1,"_",eQTL_name,"_eqtl_sum.txt --cojo-cond ",independent_eqtl_SNP_to_contition_on,"_eqtl_independent.snp --out ",variant_id,"_",qtl1,"_",eQTL_name,"_",independent_eqtl_SNP_to_contition_on,"eqtl_step2", sep=""))

                conditioned_dataset_eQTL=fread(paste0(variant_id,"_",qtl1,"_",eQTL_name,"_",independent_eqtl_SNP_to_contition_on,"eqtl_step2.cma.cojo"))
                # Conditioned dataset doesnt nontain the SNP that we condition the data to, this needs to be included.
                conditioned_dataset_eQTL_condSNP=fread(paste0(variant_id,"_",qtl1,"_",eQTL_name,"eqtl_step1.jma.cojo"))
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
