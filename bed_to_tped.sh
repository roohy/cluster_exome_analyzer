#!/usr/bin/env bash


clear
echo "Hello $USER"
now=$(date)
echo "Today is $now"

module load plink/1.90
BASE_LOCATION="/sc/arion/scratch/belbig01/for_roohy/exome/UKBexomeOQFE_chr"
BED_PRE="${BASE_LOCATION}UKBexomeOQFE_chr"
BED_SUFF="_v2.bed"



LOG_PRE="/hpc/users/shemir03/exome/plink_sub/"
SUB_FILE_PRE="SUBFILE"


for(( i = 1 ; i <= 22 ; i++))
do
sub_addr="${LOG_PRE}${SUB_FILE_PRE}_${i}"
echo "#BSUB -J ${i}_plinktped" > $sub_addr
echo "#BSUB -P acc_ipm2" >> $sub_addr
echo "#BSUB -q premium" >> $sub_addr
echo "#BSUB -n 1" >> $sub_addr
echo '#BSUB -R "span[hosts=1] affinity[core(4, same=socket, exclusive=(socket, injob))]"' >> $sub_addr
echo "#BSUB -R rusage[mem=80000]" >> $sub_addr
echo "#BSUB -R himem" >> $sub_addr
echo "#BSUB -W 4:00" >> $sub_addr
echo "#BSUB -o ${LOG_PRE}plink_chr${i}.stdout" >> $sub_addr
echo "#BSUB -e ${LOG_PRE}plink_chr${i}.stderr" >> $sub_addr
echo "#BSUB -L /bin/bash" >> $sub_addr

echo "module load plink/1.90" >> $sub_addr


echo "plink --bfile ${1}${i} --recode transpose --out ${2}_chr${i}" >> $sub_addr
bsub < $sub_addr



done