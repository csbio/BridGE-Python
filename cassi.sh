# use CASSI software to compute SNP-SNP interaction
# https://www.staff.ncl.ac.uk/richard.howey/cassi/

#!/bin/bash

# check argument
if [ "$#" -lt  "2" ]; then
	echo "gwas file or network number not provided"
fi

export CURRENTDIR=`pwd`
export PATH=$CURRENTDIR/cassi/:$PATH
chmod u+x cassi/cassi
alias cassi=cassi/cassi
#chmod u+x cassi

# find project directory
project_dir="$(dirname "${1}")"

# Step 1, compute SNP-SNP interaction using CASSI
# example CASSI command for case/control study
# logistic regression, only keep interactions with p-value <= 0.1 
infile=${1}
network=${2}

if [ "$network" == "0" ];
then
	cassi -i ${infile}.bed  -lr -lr-th 0.1 -max 0 -o ${project_dir}/ssM_cassi_LR_R0.txt
	python3 corefuns/cassissm.py ${infile}.bim  ${project_dir}/ssM_cassi_LR_R0.txt lr ${project_dir}/ssM_cassi_LR_R0.pkl
else
	plink --bfile $infile --seed 6377474 --make-perm-pheno $network --out ${project_dir}/rand_pheno$network
	plink --bfile $infile --pheno ${project_dir}/rand_pheno${network}.pphe --mpheno $network --make-bed --out ${infile}_R$network
	cassi -i ${infile}_R${network}.bed -lr -lr-th 0.1 -max 0 -o ${project_dir}/ssM_cassi_LR_R${network}.txt
	python3 corefuns/cassissm.py ${infile}.bim  ${project_dir}/ssM_cassi_LR_R${network}.txt lr ${project_dir}/ssM_cassi_LR_R${network}.pkl
	rm ${project_dir}/rand_pheno${network}.pphe
fi




