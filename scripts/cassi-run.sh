# This script uses CASSI software to compute SNP-SNP interaction and then converts it to BridGE's interaction format(Pickle file)
# Options:
#       - --plinkFile: input plink file without the extension
#       - --GI: genetic interaction type, similar to CASSI. Check CASSI for more information.
#       - --pvalue: p-value cutoff for filtering interactions in CASSI.
#       - --nSample: number of sample permutations to be used for generating random interaction networks.
#       - --id: random network number (0 means real network) to be computed. This can be used instead of --nSample for creating only one network.
# https://www.staff.ncl.ac.uk/richard.howey/cassi/

#!/bin/bash

# parsing parameters
for i in "$@"; do
        case $i in
                --plinkFile=*)
                        infile="${i#*=}"
                        shift
                        ;;
                --GI=*)
                        gi="${i#*=}"
                        shift
                        ;;
                --pvalue=*)
                        pval="${i#*=}"
                        shift
                        ;;
                --nSample=*)
                        n="${i#*=}"
                        shift
                        ;;
                --id=*)
                        net_id="${i#*=}"
                        shift
                        ;;
                --*)
                        echo "Unknown option $i"
                        exit 1
                        ;;
                *)
                        ;;
        esac
done

if [ -z "$infile" ]; then
        echo "plink input file not provided";
        exit 1;
fi

if [ -z "$gi" ]; then
        echo "CASSI interaction type not provided";
        exit 1;
fi

if [ -z "$pval" ]; then
        echo "p-value cutoff not provided";
        exit 1;
fi

if [ -z "$n" ] && [ -z "$net_id" ] ; then
        echo "number of sample permutations and network id not provided";
        exit 1;
fi



# find project directory
project_dir="$(dirname "${infile}")"


# Dciding on the command based on interaction type
case $gi in
        "lr")
                # for lr
                th_command="-lr -lr-th"
                outfile="ssM_cassi_LR"
                ;;
        "lin")
                th_command="-lin -lin-th"
                outfile="ssM_cassi_LIN"
                ;;
        "je")
                th_command="-je -je-thcc"
                outfile="ssM_cassi_JE"
                ;;
        "awu")
                th_command="-awu -awu-thcc"
                outfile="ssM_cassi_AWU"  
                ;;
        "afe")
                th_command="-afe -afe-thcc"
                outfile="ssM_cassi_AFE"
                ;;
        "wz")
                th_command="-wz -wz-thcc"
                outfile="ssM_cassi_WZ"
                ;;
        *)
                ;;
esac


seed=6377474
if [ -z "$net_id" ] ; then 

        # Step 1, compute SNP-SNP interaction using CASSI
        cassi -i ${infile}.bed  ${th_command} ${pval} -max 0 -o ${project_dir}/${outfile}_R0.txt
        # Step 2, convert SNP-SNP interaction to Python
        python3 -m cassissm ${infile}.bim  ${project_dir}/${outfile}_R0.txt ${gi} ${project_dir}/${outfile}_R0.pkl
        for i in `seq 1 $n` ;
        do
                tmp_seed=$(($seed+$i))
                plink --bfile $infile --seed $tmp_seed --make-perm-pheno 1 --out ${project_dir}/rand_pheno${i}
                plink --bfile $infile --pheno ${project_dir}/rand_pheno${i}.pphe --mpheno 1 --make-bed --out ${infile}_R${i}
                cassi -i ${infile}_R${i}.bed ${th_command} ${pval}  -max 0 -o ${project_dir}/${outfile}_R${i}.txt
                python3 -m cassissm ${infile}.bim  ${project_dir}/${outfile}_R${i}.txt ${gi} ${project_dir}/${outfile}_R${i}.pkl
                rm ${project_dir}/rand_pheno${i}.pphe
        done
else
        # find which network is being created
        if [ "$net_id" == "0" ];
        then
                cassi -i ${infile}.bed  ${th_command} ${pval} -max 0 -o ${project_dir}/${outfile}_R0.txt
                python3 -m cassissm ${infile}.bim  ${project_dir}/${outfile}_R0.txt ${gi} ${project_dir}/${outfile}_R0.pkl
        else
                tmp_seed=$(($seed+$net_id))
                plink --bfile $infile --seed $tmp_seed --make-perm-pheno 1 --out ${project_dir}/rand_pheno${net_id}
                plink --bfile $infile --pheno ${project_dir}/rand_pheno${net_id}.pphe --mpheno 1 --make-bed --out ${infile}_R${net_id}
                cassi -i ${infile}_R${net_id}.bed ${th_command} ${pval}  -max 0 -o ${project_dir}/${outfile}_R${net_id}.txt
                python3 -m cassissm ${infile}.bim  ${project_dir}/${outfile}_R${net_id}.txt ${gi} ${project_dir}/${outfile}_R${net_id}.pkl
                rm ${project_dir}/rand_pheno${net_id}.pphe
        fi

fi







