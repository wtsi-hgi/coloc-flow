#!/bin/bash


# gwasZ=($(ls *GWAS.txt))
# qtlZ=($(ls *QTL.txt))
# gwasLD=($(ls *GWAS.txt.ld))
# qtlLD=($(ls *QTL.txt.ld))
# 
# for (( i=0 ; i < ${#gwasZ[@]} ; ++i )) ;
# do
# 	echo $i
# 	echo ""
# 	echo ${gwasLD[i]}
# 	echo ${gwasZ[i]}
# 	echo ${qtlZ[i]}
# 	echo ${qtlLD[i]}
# 	/home/c/cs806/caviar/CAVIAR-C++/eCAVIAR -l ${gwasLD[i]} -z ${gwasZ[i]} -l  ${qtlLD[i]} -z ${qtlZ[i]} -o ${gwasZ[i]}
# done

for filename in *GWAS.txt; do
    fname=$(echo "${filename%GWAS.txt}")
    
    gwasZ=$(echo ${fname}GWAS.txt)
    qtlZ=$(echo ${fname}QTL.txt)
    gwasLD=$(echo ${fname}GWAS.txt.ld)
    qtlLD=$(echo ${fname}QTL.txt.ld)
    echo ""
    echo ${gwasZ}
    echo ${gwasLD}
    echo ${qtlZ}
    echo ${qtlLD}
    
    /home/c/cs806/caviar/CAVIAR-C++/eCAVIAR -l ${gwasLD} -z ${gwasZ} -l  ${qtlLD} -z ${qtlZ} -o ${gwasZ}
    
done

#/home/c/cs806/caviar/CAVIAR-C++/eCAVIAR -l ENSG00000169684_for_eCAVIAR_GWAS.txt.ld -z ENSG00000169684_for_eCAVIAR_GWAS.txt -l  ENSG00000169684_for_eCAVIAR_QTL.txt.ld -z ENSG00000169684_for_eCAVIAR_QTL.txt -o testCAVIAR2
