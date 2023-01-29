#!/bin/bash
# remove SNPs without assigned genotype

PlinkFile=$1
OutputFile=$2

awk '$5 != "0" && $6 != "0" { print $0 }' ${PlinkFile}.bim > tmp.snp
awk '$5 == "0" || $6 == "0" { print $0 }' ${PlinkFile}.bim > nogenotype_snp
DIRPATH=$(dirname "$plinkFile")
plink --bfile ${PlinkFile} --extract tmp.snp --allow-no-sex --make-bed --out ${OutputFile}

rm tmp.snp nogenotype_snp
