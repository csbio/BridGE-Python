#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: data_checkpopulation.sh <[options]>

 DESCRIPTION
 This script combines study data with 1000 project data and check the population.
 Make sure the study data uses the same build as 1000 project data.

 GENERAL INPUTS:

  --plinkFile=PLINKFILE
        PLINK filename without the extension.
        For example, your plink files are plinkexample.bim, plinkexample.bed,
        plinkexample.fam. Please use "--plinkFile=plinkexample".
        This input is required.

  --prj1000Pop=PRJ1000POP
       1000 genome project population
       There too many populations in the project
       Only CEU, CHB, YRI and ASW populations will be plotted
       The is to make sure the population of study data will be included 
       population code at https://www.internationalgenome.org/category/population/
       CEU: European
       ASW: African
       CHB: Chinese
       YRI: Youruba
       This input is optional
       Default is GIH

--prj1000File=PRJ1000FILE
       1000 genome project PLINK file
       This input is optional 
       Default is ${bridgeDir}/1000_genome_project/ALL.shapeit2_integrated_v1a.GRCh38.20181129.phased.rsid

  --popIDFile=POPIDFILE
       population ID file for 1000 project samples
       This input is optional
       Default is ${bridgeDir}/1000_genome_project/allpopid.txt
EOF
}

# get input variables
options=$@
arguments=${options}


for argument in $options
do
     case $argument in
     --plinkFile=*) plinkFile=${argument/*=/""} ;;
     --bridgeDir=*) bridgeDir=${argument/*=/""} ;;
     --legendPos=*) legendPos=${argument/*=/""} ;;
     --prj1000Pop=*) prj1000Pop=${argument/*=/""} ;;
     --prj1000File=*) prj1000File=${argument/*=/""} ;;
     --popIDFile=*) popIDFile=${argument/*=/""} ;;
     esac
done


if [ -z "${plinkFile}" ]; then printf "Please specify plinkFile."; exit; fi
if [ -z ${bridgeDir} ]; then bridgeDir=`pwd`; fi
if [ -z ${legendPos} ]; then legendPos=southeast; fi
if [ -z ${prj1000Pop} ]; then prj1000Pop=GIH; fi
if [ -z ${prj1000File} ]; then prj1000File=${bridgeDir}/1000_genome_project/ALL.shapeit2_integrated_v1a.GRCh38.20181129.phased.rsid; fi
if [ -z ${popIDFile} ]; then popIDFile=${bridgeDir}/1000_genome_project/allpopid.txt; fi



DIRPATH=$(dirname "$plinkFile")
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

# 1. merge with 1000 project data
# only keep snps from study

plink --bfile ${prj1000File} --extract raw/${plinkFile}.bim --allow-no-sex --make-bed --out intermediate/prj1000_tmp0 > /dev/null

# only keep CEU, CHB, YRI and ASW populations
grep -w CEU ${popIDFile} | awk '{print $1 "\t" $2}' > intermediate/sublist.tmp
grep -w CHB ${popIDFile} | awk '{print $1 "\t" $2}' >> intermediate/sublist.tmp
grep -w ASW ${popIDFile} | awk '{print $1 "\t" $2}' >> intermediate/sublist.tmp
grep -w YRI ${popIDFile} | awk '{print $1 "\t" $2}' >> intermediate/sublist.tmp

# add ${prj1000Pop} if it's not in CEU, CHB, YRI and ASW 
if [[ ! "${prj1000Pop}" =~ ^(CEU|CHB|ASW|YRI)$ ]]; then 
     grep -w ${prj1000Pop} ${popIDFile} | awk '{print $1 "\t" $2}' >> intermediate/sublist.tmp
fi

# filter out nonfounders
plink --bfile intermediate/prj1000_tmp0 -filter-founders --allow-no-sex --make-bed --out intermediate/prj1000_tmp1 > /dev/null

# only keep subjects in sublist.tmp
plink --bfile intermediate/prj1000_tmp1 --keep intermediate/sublist.tmp --allow-no-sex --make-bed --out intermediate/prj1000_tmp2 > /dev/null

# only keep biallelic SNPs
plink --bfile intermediate/prj1000_tmp2 --biallelic-only strict --allow-no-sex --make-bed --out intermediate/prj1000_tmp3 > /dev/null

# merge 1000 project data with study data
plink --bfile intermediate/prj1000_tmp3 --bmerge raw/${plinkFile} --allow-no-sex --make-bed --out intermediate/allpop_tmp > /dev/null

# if any SNPs in two datasets have different alleles, there will be error "variants with 3+ alleles present"
if [ -f "allpop_tmp-merge.missnp" ]; then
  plink --bfile intermediate/prj1000_tmp3 --exclude intermediate/allpop_tmp-merge.missnp --allow-no-sex --make-bed --out intermediate/prj1000_tmp4 > /dev/null
  plink --bfile raw/${plinkFile} --exclude intermediate/allpop_tmp-merge.missnp --allow-no-sex --make-bed --out  intermediate/${plinkFile}_tmp1 > /dev/null
  plink --bfile intermediate/prj1000_tmp4 --bmerge intermediate/${plinkFile}_tmp1 --allow-no-sex --make-bed --out intermediate/allpop_tmp > /dev/null
fi


# merge could fail due to mismatched SNPs
# remove mismatched SNPs and merge data again
while [ ! -f intermediate/allpop_tmp.bim ] && [ -f intermediate/allpop_tmp-merge.missnp ];
do
  plink --bfile intermediate/prj1000_tmp2 --exclude intermediate/allpop_tmp-merge.missnp --allow-no-sex --make-bed --out intermediate/plink_tmp > /dev/null
  plink --bfile raw/${plinkFile} --exclude intermediate/allpop_tmp-merge.missnp --allow-no-sex --make-bed --out intermediate/${plinkFile}_tmp  > /dev/null
  rm intermediate/allpop_tmp-merge.missnp
  plink --bfile intermediate/plink_tmp --bmerge intermediate/${plinkFile}_tmp --allow-no-sex --make-bed --out intermediate/allpop_tmp > /dev/null
done

# genereate population idex
z=0
for i in `echo CEU CHB ASW YRI ${prj1000Pop} | xargs -n1 | sort -u | xargs`; do z=`echo ${z} 0`; done

echo `echo CEU CHB ASW YRI ${prj1000Pop} | xargs -n1 | sort -u | xargs` StudyPop > intermediate/allpopid.txt

while read line
do
  pattern=`echo ${line} | awk '{print $1 " " $2}'`
  pop=`grep -w "${pattern}" ${popIDFile} | awk '{print $3}'`

  if [ "${pop}" == "CEU" ]; then
    echo ${pattern} `echo $z | awk '{ $1=1; print }'` >> intermediate/allpopid.txt
  elif [ "${pop}" == "CHB" ]; then
    echo ${pattern} `echo $z | awk '{ $2=1; print }'` >> intermediate/allpopid.txt
  elif [ "${pop}" == "ASW" ]; then
    echo ${pattern} `echo $z | awk '{ $3=1; print }'` >> intermediate/allpopid.txt
  elif [ "${pop}" == "${prj1000Pop}" ]; then
    echo ${pattern} `echo $z | awk '{ $4=1; print }'` >> intermediate/allpopid.txt
  elif [ "${pop}" == "YRI" ]; then
    echo ${pattern} `echo $z | awk '{ $5=1;print }'` >> intermediate/allpopid.txt
  else
    echo ${pattern} `echo $z | awk -v x=6, '{ $x=1; print }'` >> intermediate/allpopid.txt
  fi
done < intermediate/allpop_tmp.fam

# To perform principal components or MDS analysis, it is also important that we use SNPs that are not too correlated.
plink --bfile intermediate/allpop_tmp --allow-no-sex --indep-pairwise 50 5 0.2 --out intermediate/prunedsnps_tmp > /dev/null
plink --bfile intermediate/allpop_tmp --allow-no-sex --extract intermediate/prunedsnps_tmp.prune.in --make-bed --out intermediate/prunedsnps_tmp > /dev/null

# Calculate genome-wide estimates of IBD sharing
plink --bfile intermediate/prunedsnps_tmp --freq --allow-no-sex --out intermediate/prunedsnps_tmp > /dev/null
plink --bfile intermediate/prunedsnps_tmp --read-freq intermediate/prunedsnps_tmp.frq --genome --allow-no-sex --out intermediate/prunedsnps_tmp > /dev/null

# Compute MDS scores
plink --bfile intermediate/prunedsnps_tmp --allow-no-sex --read-genome intermediate/prunedsnps_tmp.genome --cluster --mds-plot 2 --out intermediate/prunedsnps_tmp > /dev/null

cp intermediate/prunedsnps_tmp.mds intermediate/${plinkFile}_prj1000.mds
python3 -m plotmds intermediate/${plinkFile}_prj1000.mds intermediate/allpopid.txt ${prj1000Pop} intermediate/MDS_${plinkFile}

rm intermediate/*tmp*
rm intermediate/allpopid.txt





