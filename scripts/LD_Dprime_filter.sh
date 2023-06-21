#!/bin/bash

# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --plinkFile=*) plinkFile=${argument/*=/""} ;;
        --ldWindow=*) ldWindow=${argument/*=/""} ;;
        --ldDprime=*) ldDprime=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

if [ -z "${plinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${ldDprime}" ]; then ldDprime=0.5; fi
if [ -z "${ldWindow}" ]; then ldWindow=50; fi


# filter based on Dprime if dprime value is set
# Compute Dprime values
plink --bfile ${plinkFile} --r2 dprime with-freqs --ld-window-r2 0.001 --ld-window ${ldWindow}
# choose SNP with lower MAF from pairs with D' > ldDprime to be removed
awk -v D=${ldDprime} 'NR > 1 && $10 > D { if ($4 < $8) {print $3} else {print $7} }' plink.ld > rm_plink
# remove the list from the data
plink --bfile ${plinkFile} --exclude rm_plink --make-bed --out ${plinkFile}_dprime
plink --bfile ${plinkFile}_dprime --allow-no-sex --noweb --recodeA --out ${plinkFile}_dprime
rm rm_plink
rm plink.*
