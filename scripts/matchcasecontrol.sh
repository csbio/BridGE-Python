

#!/bin/bash

# this script is used to match cases and controls using plink's size-2 clustering
plinkFile=$1  # input plink file
outputFile=$2 # new plink file with matched cases and controls 

# compute ibd file
if [ ! -f ${plinkFile}.genome ];
then
        ./plink --bfile ${plinkFile} --freq --allow-no-sex --out ${plinkFile}
        ./plink --bfile ${plinkFile} --read-freq ${plinkFile}.frq --genome --allow-no-sex --out ${plinkFile}
fi

# size2 clustering
./plink --bfile ${plinkFile} --read-genome ${plinkFile}.genome --cluster --cc --mc 2 --ppc 0.01 --noweb --allow-no-sex --out ${plinkFile}

# get subject list based on clustering
cat  ${plinkFile}.cluster2|awk '{print $3}'| sort  | uniq -c | sort -nr > cluster2_count.tmp
awk '$1==1' cluster2_count.tmp|awk '{print $2}' > sol.tmp
cp ${plinkFile}.cluster2 cluster2.tmp
for val in `cat sol.tmp`
do
	awk -v x=$val '$3!=x' cluster2.tmp > tmp.tmp
	mv tmp.tmp cluster2.tmp 
done
mv ${plinkFile}.cluster2 ${plinkFile}.cluster2.orig
mv cluster2.tmp  ${plinkFile}.cluster2 
cat ${plinkFile}.cluster2  | awk '{print $1 "\t" $2 }' >  cluster2_subject2keep

# generate new plink data
./plink --bfile ${plinkFile} --allow-no-sex --keep cluster2_subject2keep --make-bed --out ${outputFile}

rm cluster2_count.tmp sol.tmp
