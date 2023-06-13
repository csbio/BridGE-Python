#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: data_remove_outlier.sh <[options]>

 DESCRIPTION
 Data processing version 1 (V1_NIND) : no imputation, no deleterious SNP prioritization

 GENERAL INPUTS:

  --plinkFile=PLINKFILE
        PLINK filename without the extension.
        For example, your plink files are plinkexample.bim, plinkexample.bed,
        plinkexample.fam. Please use "--plinkFile=plinkexample".
        This input is required.
 
  --mdsFile=MDSFILE
        mds file genereated from PLINK
        This input is required

  --x1=X1
        Samples need to have C1 greater than cutoffX1
        Samples don't pass this criteria will be filtered out
        Ths input is required

  --x2=X2
        Samples need to have C1 smaller than cutoffX1
        Samples don't pass this criteria will be filtered out
        This input is required

  --y1=Y1
        Samples need to have C2 greater than cutoffY1
        Samples don't pass this criteria will be filtered out 
        This input is required

  --y2=Y2
        Samples need to have C2 smaller than cutoffY1
        Samples don't pass this criteria will be filtered out
        This input is required
EOF
}

# get input variables
options=$@
arguments=${options}

for argument in $options
do
     case $argument in
     --plinkFile=*) plinkFile=${argument/*=/""} ;;
     --mdsFile=*) mdsFile=${argument/*=/""} ;;
     --x1=*) x1=${argument/*=/""} ;;
     --x2=*) x2=${argument/*=/""} ;;
     --y1=*) y1=${argument/*=/""} ;;
     --y2=*) y2=${argument/*=/""} ;;
     esac
done

# check variables
if [ -z "${plinkFile}" ]; then printf "Please specify plinkFile."; exit; fi
if [ -z "${mdsFile}" ]; then printf "Please specify MDS file"; exit; fi
if [ -z "${x1}" ] || [ -z "${x2}" ] || [ -z "${y1}" ] || [ -z "${y2}" ]; then 
     printf "Missing inputs for x1, x2, y1, or y2, please check usage"; exit; 
fi

DIRPATH=$(dirname "$plinkFile")
plinkFile=${plinkFile##*/}
mdsFile=${mdsFile##*/}
cd $DIRPATH

awk -v X1=${x1} -v X2=${x2} -v Y1=${y1} -v Y2=${y2} '$4 >= X1 && $4 <= X2 && $5 >= Y1 && $5 <= Y2 ' ${mdsFile} | awk '{print $1,$2}' > tmp_sample2keep

plink --bfile ${plinkFile} --keep tmp_sample2keep --allow-no-sex --make-bed --out ${plinkFile}.rmoutlier

rm tmp_sample2keep




