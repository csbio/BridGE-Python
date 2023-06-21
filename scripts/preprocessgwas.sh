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
#if [ -z "${ldDprime}" ]; then ldDprime=0.1; fi



DIRPATH=$(dirname "$plinkFile")
# now using raw/ sub directory
DIRPATH=$(dirname "$DIRPATH")
plinkFile=${plinkFile##*/}
cd $DIRPATH

# create sub directories if do not alredy exist

if [ ! -d "intermediate" ]; then
     mkdir intermediate
fi

if [ ! -d "results" ]; then
     mkdir results
fi


plink --bfile raw/${plinkFile} --noweb --allow-no-sex --mind ${mind} --geno ${geno} --maf ${maf} --hwe ${hwe} --make-bed --out intermediate/${plinkFile}_tmp0

extractchr1-22.sh intermediate/${plinkFile}_tmp0 intermediate/${plinkFile}_tmp1

excludenogenotypesnps.sh intermediate/${plinkFile}_tmp1 intermediate/${plinkFile}_tmp2
removerelatedindividual.sh intermediate/${plinkFile}_tmp2 intermediate/${plinkFile}_tmp3 ${pihat}
if [ "${matchCC}" -eq 1 ] 
then
     matchcasecontrol.sh intermediate/${plinkFile}_tmp3 intermediate/gwas_data_all 
     mv intermediate/${plinkFile}_tmp3.cluster1 intermediate/PlinkFile.cluster1
     mv intermediate/${plinkFile}_tmp3.cluster2 intermediate/PlinkFile.cluster2
     mv intermediate/${plinkFile}_tmp3.cluster2.orig intermediate/PlinkFile.cluster2.orig
else
     plink --bfile intermediate/${plinkFile}_tmp3 --make-bed --out intermediate/gwas_data_all
fi



# get less redundant SNP set based on R2
plink --bfile intermediate/gwas_data_all --allow-no-sex \
     --indep-pairwise ${ldWindow} ${ldShift} ${ldR2} --noweb \
     --out intermediate/gwas_data_all

# generate new plink data
plink --bfile intermediate/gwas_data_all --allow-no-sex --extract intermediate/gwas_data_all.prune.in \
     --make-bed --out intermediate/gwas_data_final


# convert plink to python readable
plink --bfile intermediate/gwas_data_final --allow-no-sex --noweb --recodeA \
	--out intermediate/RecodeA_file

# Compute all R2 pirwise data to later be used in get_interaction_list 
plink --bfile intermediate/gwas_data_final --r2 square

mv intermediate/RecodeA_file.raw intermediate/gwas_data_final.raw
mv plink.ld intermediate/plink.ld



rm intermediate/*tmp*
rm intermediate/gwas_data_all.*
rm plink*
