### Workflow

##Extract benchmarking data if not do this before
#grep ">" all_1224.fasta > all_1224_1.fasta
#sed 's/.//' all_1224_1.fasta > all_1224_TRM.txt
#cat all_1224_TRM.txt | sort | uniq | sort > all_1224.txt
#rm all_1224_1.fasta
#rm all_1224_TRM.txt


Iteration=1 # Setting iteration number 
JackhmmerEvalue=15 # Setting Evalue as jackhmmer threshold
CutoffEvalue=30   # Setting Evalue as cutoff

## Processing progress
for f in $(ls *.fasta); do           # Please ensure every query sequences are protein sequences and in the form of "fasta"
echo "Processing file: ${f}"; 
gene=`echo ${f%??????}`;

##Simulating gene family by uniref50
jackhmmer -N ${Iteration} -A ${gene}.sto -E 1e-${JackhmmerEvalue}  --cpu 5 --noali --notextw --acc --chkhmm ${gene}  ${gene}.fasta uniref50.fas > ${gene}_summary.txt  # Please ensure references database in the form of "fas"


##Construct HMM 
#Change the number when set different iteration number
#hmmbuild ${gene}-${Iteration}.hmm ${gene}.sto
hmmpress ${gene}-${Iteration}.hmm

##Benchmarking searching
hmmsearch --tblout ${gene}_hmm.txt -E 1e-${CutoffEvalue} --cpu 5 --noali --notextw --acc ${gene}-${Iteration}.hmm all_1224.fas > ${gene}_hmm_summary.txt # Please ensure benchmarking database in the form of "fas"

##Trimming Searching result

sed -i '1,3d' ${gene}_hmm.txt 
A=$(sed -n '$=' ${gene}_hmm.txt)
sed -i $(($A-10+1)),${A}d ${gene}_hmm.txt 
awk '{print $1}' ${gene}_hmm.txt > ${gene}_hmm_T.txt

cat ${gene}_hmm_T.txt | sort | uniq | sort > ${gene}_hmm_T_u.txt  
diff  all_1224.txt  ${gene}_hmm_T_u.txt  | grep \<  | awk ' $1 = " " ' > ${gene}_hmm_F.txt 

##Extract the TP,TN,FP,FN
grep ${gene} ${gene}_hmm_T_u.txt > ${gene}TP.txt
cat ${gene}TP.txt | sort | uniq | sort > ${gene}TP_u.txt  
diff  ${gene}_hmm_T_u.txt  ${gene}TP_u.txt  | grep \<  | awk ' $1 = " " ' > ${gene}FP.txt 

grep ${gene} ${gene}_hmm_F.txt > ${gene}FN.txt
cat ${gene}FN.txt | sort | uniq | sort > ${gene}FN_u.txt
diff  ${gene}_hmm_F.txt  ${gene}FN_u.txt  | grep \<  | awk ' $1 = " " ' > ${gene}TN.txt 


##Count the line which represents the number of sample

touch ${gene}_Result.txt

TP=`awk 'END{print NR}' ${gene}TP.txt` 
TN=`awk 'END{print NR}' ${gene}TN.txt`
FP=`awk 'END{print NR}' ${gene}FP.txt` 
FN=`awk 'END{print NR}' ${gene}FN.txt`
echo "TP:${TP}" >> ${gene}_Result.txt
echo "TN:${TN}" >> ${gene}_Result.txt
echo "FP:${FP}" >> ${gene}_Result.txt
echo "FN:${FN}" >> ${gene}_Result.txt

##Remove the file exluding results
#Want to see the process just need to add '#' before every lines
rm ${gene}.sto
rm ${gene}_hmm.txt
rm ${gene}_hmm_F.txt
rm ${gene}_hmm_summary.txt
rm ${gene}_hmm_T.txt
rm ${gene}_hmm_T_u.txt
rm ${gene}_summary.txt
rm ${gene}-1.hmm
rm ${gene}-1.hmm.h3f
rm ${gene}-1.hmm.h3i
rm ${gene}-1.hmm.h3m
rm ${gene}-1.hmm.h3p
rm ${gene}-2.hmm
rm ${gene}-2.hmm.h3f
rm ${gene}-2.hmm.h3i
rm ${gene}-2.hmm.h3m
rm ${gene}-2.hmm.h3p
rm ${gene}-3.hmm
rm ${gene}-3.hmm.h3f
rm ${gene}-3.hmm.h3i
rm ${gene}-3.hmm.h3m
rm ${gene}-3.hmm.h3p
rm ${gene}FN.txt
rm ${gene}FN_u.txt
rm ${gene}FP.txt
rm ${gene}TN.txt
rm ${gene}TP.txt
rm ${gene}TP_u.txt

echo "${f} file has been processed, mission complete";
done

##End