This files required in input are:
1) GWAS summary statistics that contain these fields (colnames may be present): "rsID", "CHR", "BP", "NEA", "EA", "EAF", "BETA", "SE", "Z", "P", "N", "CASES", "CONTROLS", "PHENO"
2) Full eQTL summary statistics that contain these fields in the right order with no header: "gene", "SNP", "distance_to_TSS", "p", "beta"
The SNPs have to contain rsids. 
3) plink file .fam .bed .bim that ideally has variants annotated with rs'ids

