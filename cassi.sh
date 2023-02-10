# use CASSI software to compute SNP-SNP interaction
# https://www.staff.ncl.ac.uk/richard.howey/cassi/

#!/bin/bash

# check argument
if [ "$#" -lt  "1" ]; then
        echo "gwas file not provided"
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
infile="${1}.bed"
cassi -i $infile  -lr -lr-th 0.1 -max 0 -o ${project_dir}/ssM_cassi_LR_R0.txt

# joint effect, only keep interactions with p-value <= 0.1
cassi -i $infile -je -je-thcc 0.1 -max 0 -o ${project_dir}/ssM_cassi_JE_R0.txt

# adjusted Wu method, only keep interactions with p-value <= 0.1
cassi -i $infile -awu -awu-thcc 0.1 -max 0 -o ${project_dir}/ssM_cassi_AWU_R0.txt

# Adjusted Fast Epistasis, only keep interactions with p-value <= 0.1
cassi -i $infile -afe -afe-thcc 0.1 -max 0 -o ${project_dir}/ssM_cassi_AFE_R0.tx

# Wellek Ziegler method, only keep interactions with p-value <= 0.1
cassi -i $infile -wz -wz-thcc 0.1 -max 0 -o ${project_dir}/ssM_cassi_WZ_R0.txt

# for quantitative traits
cassi -i $infile -lin -lin-th 0.1 -max 0 -o ${project_dir}/ssM_cassi_LIN_R0.txt

python3 corefuns/cassissm.py $infile  ${project_dir}/ssM_cassi_LR_R0.txt lr ${project_dir}/ssM_cassi_LR_R0.pkl

plink --bfile $infile --make-perm-pheno 10 --out ${project_dir}/rand_pheno

for R in `seq 1 10`
do
  plink --bfile $infile --pheno ${project_dir}/rand_pheno.pphe --mpheno $R --make-bed --out ${infile}_R$R
done

cassi -i ${infile}_R1.bed -lr -lr-th 0.1 -max 0 -o ${project_dir}/ssM_cassi_LR_R1.txt
python3 corefuns/cassissm.py $infile  ${project_dir}/ssM_cassi_LR_R1.txt lr ${project_dir}/ssM_cassi_LR_R1.pkl