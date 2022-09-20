#!/bin/bash

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

./plink --bfile ${plinkFile} --noweb --allow-no-sex --mind ${mind} --geno ${geno} --maf ${maf} --hwe ${hwe} --make-bed --out ${plinkFile}_tmp0

scripts/extractchr1-22.sh ${plinkFile}_tmp0 ${plinkFile}_tmp1

scripts/excludenogenotypesnps.sh ${plinkFile}_tmp1 ${plinkFile}_tmp2
scripts/removerelatedindividual.sh ${plinkFile}_tmp2 ${plinkFile}_tmp3 ${pihat}
if [ "${matchCC}" -eq 1 ] 
then
     scripts/matchcasecontrol.sh ${plinkFile}_tmp3 gwas_data_all 
     mv ${plinkFile}_tmp3.cluster1 PlinkFile.cluster1
     mv ${plinkFile}_tmp3.cluster2 PlinkFile.cluster2
     mv ${plinkFile}_tmp3.cluster2.orig PlinkFile.cluster2.orig
else
     ./plink --bfile ${plinkFile}_tmp3 --make-bed --out gwas_data_all
fi


# get less redundant SNP set
./plink --bfile gwas_data_all --allow-no-sex \
     --indep-pairwise ${ldWindow} ${ldShift} ${ldR2} --noweb \
     --out gwas_data_all

# generate new plink data
./plink --bfile gwas_data_all --allow-no-sex --extract gwas_data_all.prune.in \
	--make-bed --out gwas_data_final


# conver plink to python readable
./plink --bfile gwas_data_final --allow-no-sex --noweb --recodeA \
	--out RecodeA_file

mv gwas_data_final data/
mv data/RecodeA_file.raw data/gwas_data_final.raw
