#!/bin/sh
module load HMMER/3.3.2-gompi-2021a;
####
##number of iters
n=1
#p 
p=4
for f in $(ls Genes/*.fasta); do
echo "Processing file: ${f}";
jackhmmer -N $n -A ${f}.sto --incE 0.0001 --cpu $p ${f} SourceFiles/VFDB_setB_pro.fas;
done
mv Genes/*.sto ../AlnFiles/*
####
for f in $(ls AlnFiles/*); do
echo "building HMM: ${f}";
hmmbuild Profile_HMM/${f}_profile.hmm ${f};
done
cat Profile_HMM/*_profile.hmm > Yersinia_HMM.dat
hmmpress Yersinia_HMM.dat
