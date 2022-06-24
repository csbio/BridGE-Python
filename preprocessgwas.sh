#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
Usage: processgwas.sh <[options]>
DESCRIPTION  
This script is used to process data (from plink format) and prepare for 
running BridGE algorithm.

REQUIRED OPTIONS
  --plinkFile=PLINKFILE
	PLINK filename without the extension. 
	For example, your plink files are plinkexample.bim, plinkexample.bed, 
	plinkexample.fam. Please use "--plinkFile=plinkexample".

OPTIONAL OPTIONS
  --mind=MIND
	PLINK parameter used to control missing genotypes per individual. 
	Default is 0.02.

  --geno=GENO
	PLINK parameter used to control missing genotypes per SNP. 
	Default is 0.02. 

  --maf=MAF 
	PLINK parameter used to control allele frequency. 
	Default is 0.05.

  --hwe=HWE
	PLINK parameter used to control Hardy-Weinberg equilibrium. 
	Default is 0.000001.

  --pihat=PIHAT
	PIHAT: Identical twins, and duplicates, are 100% identical by descent (Pihat 1.0); 
		First-degree relatives are 50% IBD (Pihat 0.5); 
		Second-degree relatives are 25% IBD (Pihat 0.25); 
		Third-degree relatives are 12.5% equal IBD (Pihat 0.125); 
		4th degree Pihat = 0.0625; 5th degree Pihat=0.03125 level.        
	Default is 0.05.
	
  --matchCC=MATCHCC
	This parameter tells if cases and controls need to be matched.
	1: matched; 0: unmatched. 
	Default is 1.

  --ldWindow=LDWINDOW
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 50.

  --ldShift=LDSHIFT
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 5.

  --ldR2=LDR2
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 0.1. 

EOF
}


# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --plinkFile=*) plinkFile=${argument/*=/""} ;;
        --geno=*) geno=${argument/*=/""} ;;
        --mind=*) mind=${argument/*=/""} ;;
        --maf=*) maf=${argument/*=/""} ;;
        --hwe=*) hwe=${argument/*=/""} ;;
        --pihat=*) pihat=${argument/*=/""} ;;
        --matchCC=*) matchCC=${argument/*=/""} ;;
        --ldWindow=*) ldWindow=${argument/*=/""} ;;
        --ldShift=*) ldShift=${argument/*=/""} ;;
        --ldR2=*) ldR2=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${plinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${mind}" ]; then mind=0.02; fi
if [ -z "${geno}" ]; then geno=0.02; fi
if [ -z "${maf}" ]; then maf=0.05; fi
if [ -z "${hwe}" ]; then hwe=0.000001; fi
if [ -z "${pihat}" ]; then pihat=0.2; fi
if [ -z "${matchCC}" ]; then matchCC=1; fi
if [ -z "${ldWindow}" ]; then ldWindow=50; fi
if [ -z "${ldShift}" ]; then ldShift=5; fi
if [ -z "${ldR2}" ]; then ldR2=0.1; fi

# qc based on mind geno maf and hwe
plink1.9 --bfile ${plinkFile} --noweb --allow-no-sex --mind ${mind} \
	--geno ${geno} --maf ${maf} --hwe ${hwe} --make-bed --out ${plinkFile}_tmp0 > /dev/null

# only keep autosomal chromosomes
extractchr1-22.sh ${plinkFile}_tmp0 ${plinkFile}_tmp1

# remove snps without correct genotype assigned
excludenogenotypesnps.sh ${plinkFile}_tmp1 ${plinkFile}_tmp2

# remove related individuals
removerelatedindividual.sh ${plinkFile}_tmp2 ${plinkFile}_tmp3 ${pihat}

# use size 2 clustering to match case and control
if [ "${matchCC}" -eq 1 ] 
then
     matchcasecontrol.sh ${plinkFile}_tmp3 gwas_data_all 
     mv ${plinkFile}_tmp3.cluster1 PlinkFile.cluster1
     mv ${plinkFile}_tmp3.cluster2 PlinkFile.cluster2
     mv ${plinkFile}_tmp3.cluster2.orig PlinkFile.cluster2.orig
else
     plink1.9 --bfile ${plinkFile}_tmp3 --make-bed --out gwas_data_all
fi

# get less redundant SNP set
plink1.9 --bfile gwas_data_all --allow-no-sex \
     --indep-pairwise ${ldWindow} ${ldShift} ${ldR2} --noweb \
     --out gwas_data_all > /dev/null

# generate new plink data
plink1.9 --bfile gwas_data_all --allow-no-sex --extract gwas_data_all.prune.in \
	--make-bed --out gwas_data_final  > /dev/null

# conver plink to matlab
plink1.9 --bfile gwas_data_final --allow-no-sex --noweb --recodeA \
	--out RecodeA_file  > /dev/null

mv RecodeA_file.raw data/gwas_data_final.raw
mv gwas_data_final.* data/
rm *tmp*

