# This script takes three arguments
# - the path to the clean gwas summary statistics file
# - path to the clean eQTL data
# - prefix name of output file
# This script outputs a hg38 version of the gwas summary statistics file 


# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(qqman))


# List options
option_list = list(
  make_option(c("-e", "--qtl"), type="character", default=NULL, 
              help="path to prepped qtl summary statistics", metavar="character"),
  make_option(c("-g", "--gwas"), type="character", default=NULL, 
              help="path to clean gwas summary statistics", metavar="character"),
  make_option(c("-a", "--appel"), type="character", default=NULL, 
              help="prefix for output file. eg vsmcImputeEQTL_HowsonCAD2017", metavar="character"),
  make_option(c("-p", "--pcutoffQTL"), type="numeric", default = 5e-08, 
              help="pvalue threshold for selecting associations.", metavar="numeric"),
  make_option(c("-q", "--pcutoffGWAS"), type="numeric", default = 5e-05, 
              help="pvalue threshold for selecting associations.", metavar="numeric"),
  make_option(c("-c", "--candGenes"), type="character", default=NULL, 
              help="path to file containing ensembl ID of cand genes 1 per line", metavar="character")
); 


opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);

# Data input error messages
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Please provide path to summ stats of qtl.", call.=FALSE)
}

if (is.null(opt$gwas)){
  print_help(opt_parser)
  stop("Please provide path to summ stats of gwas.", call.=FALSE)
}

if (is.null(opt$appel)){
  print_help(opt_parser)
  stop("Please provide prefix name for output files.", call.=FALSE)
}


# Assign input variables --------------------------------------------------
gwasFileName <- paste0(opt$gwas)
#gwasFileName <- "../colocalization/cleanGWAS_Summary_Stats/GWAS_Alanine_Aminotransferase_Levels_Sinnott_ALT_2021_NatGenet_hg38.txt"
qtlFileName <- paste0(opt$qtl)
#qtlFileName <- "../pubAvail_QTL/GTEx/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs_pval005.txt.gz"
outFileName <- paste0(opt$appel)
#outFileName <- "vsmcImputeEQTL_HowsonCAD2017"
pcut <- as.numeric(opt$pcutoffQTL)
qcut <- as.numeric(opt$pcutoffGWAS)

cat(" \n")
cat("Started processing ", outFileName, "...... \n")
cat(" \n")

# import summary stats ------------------------------------
gwas <- data.frame(fread(gwasFileName))
gwas <- setnames(gwas, c("permID", "CHR"), c("snpid", "chr"), skip_absent=TRUE)
qtl <- data.frame(fread(qtlFileName))
qtl <- setnames(qtl, "pos", "position", skip_absent=TRUE)
qtl$snps <- sub("chr", "", qtl$snps)
invisible(gc())

# Select snps with pvalue < pcut
qtl <- qtl[which(qtl$pvalue < pcut), ]
gwas <- gwas[which(gwas$pvalue < qcut), ]
invisible(gc())

folder <- paste0("eCAVIAR_Results/eCAVIAR_", outFileName, "/")

if (dir.exists(folder)) {
  unlink(folder, recursive = TRUE)
  dir.create(folder)
}else{
  dir.create(folder)
}

# Subset genes that have genome-wide significant snp associations
sigQTL <- qtl[which(qtl$pvalue < 5e-08), ]
sigQTL <- sigQTL[order(sigQTL$pvalue), ]

# How many genes have genome wide significant associations
topSigQTL <- sigQTL[!duplicated(sigQTL$gene), ]

if (!is.null(opt$candGenes)){
  candFileName <- paste0(opt$candGenes)
  cand <- readLines(paste0(candFileName))
  topSigQTL <- qtl[qtl$gene %in% cand, ]
  topSigQTL <- topSigQTL[!duplicated(topSigQTL$gene), ]
}

# Prepare genes with genomewide significant QTL for eCAVIAR 
foreach(rr = 1:nrow(topSigQTL), .combine = rbind) %dopar% {
  eGene <- topSigQTL[rr, "gene"]
  topSNP <- topSigQTL[rr, "snps"]
  topSNP_bp <- as.numeric(topSigQTL[rr, "position"])
  topSNP_chr <- topSigQTL[rr, "chr"]
  
  # Get Summary Stats for gene
  eGene_SS <- subset(qtl, qtl$gene == eGene)
  getWin <- subset(eGene_SS, eGene_SS$position >= topSNP_bp-500000 & eGene_SS$position <= topSNP_bp+500000)
  getWin <- getWin[, c("snps", "statistic")]
  
  # Extract and match gene summ stats to gwas summ stats
  gwasMatch <- subset(gwas, gwas$hg38_markername %in% getWin$snps)
  gwasMatch <- gwasMatch[, c("hg38_markername", "zScore")]
  
  getWin <- getWin[!duplicated(getWin$snps), ]
  gwasMatch <- gwasMatch[!duplicated(gwasMatch$hg38_markername), ]
  
  getWin <- subset(getWin, getWin$snps %in% gwasMatch$hg38_markername)
  gwasMatch <- subset(gwasMatch, gwasMatch$hg38_markername %in% getWin$snps)
  
  if (nrow(getWin) < 5) {
    cat(rr, eGene, " has ", nrow(getWin), "SNPs and not eligible for eCAVIAR ...\n" )
  } else if (nrow(getWin) >= 5) {
    
    
    gwasRegList <- data.frame(gwasMatch[, "hg38_markername"])
    snpListFile <- paste0("check_in_1kG_snpList.txt")
    write.table(gwasRegList, file = snpListFile, quote = F, row.names = F, col.names = F)
    axe1 <- system(paste0("grep -f ", snpListFile, " /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge.bim  | awk \'{print $2}\'"), intern = T)
    gwasMatch <- gwasMatch[gwasMatch$hg38_markername %in% axe1, ]
    getWin <- subset(getWin, getWin$snps %in% gwasMatch$hg38_markername)
    
    # Order both dataframes so the snps are aligned
    getWin <- getWin[order(getWin$snps), ]
    gwasMatch <- gwasMatch[order(gwasMatch$hg38_markername), ]
    
    cat("Writing ", nrow(getWin), "SNPs for ", eGene, "QTL eCAVIAR analysis ...\n")
    
    qtlInput <- paste0(folder, eGene, "_for_eCAVIAR_QTL.txt")
    gwasInput <- paste0(folder, eGene, "_for_eCAVIAR_GWAS.txt")
    
    write.table(getWin, file = qtlInput, quote = F, col.names = F, row.names = F)
    write.table(gwasMatch, file = gwasInput, quote = F, col.names = F, row.names = F)
    system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge --extract ", gwasInput, " --r --matrix --threads 6 --out ", gwasInput))
    system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge --extract ", qtlInput, " --r --matrix --threads 6 --out ", qtlInput))
    system(paste0("rm ", folder,  "*.log"))
    system(paste0("rm ", folder,  "*.nosex"))
    
    
  }
  
  return(NULL)  
}


# Run eCAVIAR -------------------------------------------------------------
# system(paste0("cp /scratch/vasccell/cs806/colocalization/colocalization_scripts/eCAVIAR_script.sh ", folder))
setwd(folder)
system(paste0("bash eCAVIAR_script.sh"))
setwd("../../")


# Import and process eCAVIAR results --------------------------------------
clpp <- list.files(folder, "_col", full.names = T)

if (length(clpp) != 0) {
  clppFiles <- lapply(clpp, read.delim)
  clpp <- sub(paste0(folder, "/"), "", clpp)
  clpp <- sub("_for_eCAVIAR_GWAS.txt_col", "", clpp)
  names(clppFiles) <- clpp
  
  clppDF <- do.call(rbind, clppFiles)
  clppDF$GeneID <- sub("\\..*", "", rownames(clppDF))
  clppDF <- merge(clppDF, gwas[, c("hg38_markername", "hg38_bp", "chr", "snpid")], by.x = "SNP_ID", by.y = "hg38_markername")
  clppDF$chr <- as.numeric(clppDF$chr)
  
  #manhattan plot
  png(filename = paste0("eCAVIAR_Results/eCAVIAR_Manhattans/", outFileName, "_eCAVIARmanhattan.png"))
  
  manhattan(clppDF, chr = "chr", bp = "hg38_bp", p = "CLPP", snp = "snpid",
            col = c("gray10", "gray60"), chrlabs = NULL, logp = FALSE,
            suggestiveline = 0.05, genomewideline = FALSE,
            annotateTop = TRUE, annotatePval = 0.5, highlight = NULL)
  dev.off()
  # dev.copy(png, file = paste0("eCAVIAR_Results/eCAVIAR_Manhattans/", outFileName, "_eCAVIARmanhattan.png"))
  # dev.off()
  # dev.copy(pdf, file = paste0("eCAVIAR_Results/eCAVIAR_Manhattans/", outFileName, "_eCAVIARmanhattan.pdf"))
  # dev.off()
  
  
  sigCLPP <- subset(clppDF, CLPP >= 0.05)
  
  # Add gene names
  annot <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Annot.txt"))
  
  clppDF <- merge(clppDF, annot, by = "GeneID")
  clppDF <- clppDF[, c(1,8,7,2,3,4)]
  
  sigCLPP <- merge(sigCLPP, annot, by = "GeneID")
  sigCLPP <- sigCLPP[, c(1,8,7,2,3,4)]
  
  
  write.csv(sigCLPP, file = paste0("eCAVIAR_Results/eCAVIAR_ResTables/sig_eCAVIAR_", outFileName, ".csv"), row.names = F)
  write.csv(clppDF, file = paste0("eCAVIAR_Results/eCAVIAR_ResTables/all_eCAVIAR_", outFileName, ".csv"), row.names = F)
  
  
  cat(" \n")
  cat("Finished processing ", outFileName, "...... \n")
  cat(" \n")
  
} else if (length(clpp) == 0) {
  cat("  ", folder, " does not have clpp files")
}


