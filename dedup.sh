#!/bin/bash

# This script can be used to identify duplicates in peaks as signal
# and generate BAM with non-duplicate reads and signal duplicates

#####################################################
#### Please provide the path to dedup.sh blow ####
#####################################################

DEDUP_TOOL=/PATH to dedup.sh/



if [[ $# -ge 1 && $# -lt 5 ]]; then
    echo ""
    echo 1>&2 Error: "the required files were not provided or parameters not defined"
    echo ""
    echo 1>&2 Usage: $0 "WORK_DIR=DIRECTORY_holding_BAM_and_peak_file" "BAM=alignment_file" "PEAK=peak_file" "MAPQ=minimum_mapping_quality_score" "maxDUP=maximum_number_of_duplicates_kept_per_position" 
    echo ""
    echo "example: ${DEDUP_TOOL}"/"dedup.sh WORK_DIR=/Volumes/Users/James/ChIPseq BAM=test.bam PEAK=test.encodePeak MAPQ=20 maxDUP=5"
    echo ""
    exit 1

elif [[ $# -eq 0 ]]; then
    echo ""
    echo 1>&2 Usage: $0 "WORK_DIR=DIRECTORY_holding_BAM_and_peak_file" "BAM=alignment_file" "PEAK=peak_file" "MAPQ=minimum_mapping_quality_score" "maxDUP=maximum_number_of_duplicates_kept_per_position" 
    echo ""
    echo "example: ${DEDUP_TOOL}"/"dedup.sh WORK_DIR=/Volumes/Users/James/ChIPseq BAM=test.bam PEAK=test.encodePeak MAPQ=20 maxDUP=5"
    echo ""
    exit 1

fi


####### parse running variable

WORK_DIR=$(echo $* |awk '{gsub(/[ ]{1,},/,",",$0); gsub(/,[ ]{1,}/,",",$0); print}' |tr -s " \+" "\n" |grep -i "^WORK_DIR=" |perl -pe 's/WORK_DIR=//i')
BAM=$(echo $* |awk '{gsub(/[ ]{1,},/,",",$0); gsub(/,[ ]{1,}/,",",$0); print}' |tr -s " \+" "\n" |grep -i "^BAM=" |perl -pe 's/BAM=//i')
PEAK=$(echo $* |awk '{gsub(/[ ]{1,},/,",",$0); gsub(/,[ ]{1,}/,",",$0); print}' |tr -s " \+" "\n" |grep -i "^PEAK=" |perl -pe 's/PEAK=//i')
MAPQ=$(echo $* |awk '{gsub(/[ ]{1,},/,",",$0); gsub(/,[ ]{1,}/,",",$0); print}' |tr -s " \+" "\n" |grep -i "^MAPQ=" |perl -pe 's/MAPQ=//i')
maxDUP=$(echo $* |awk '{gsub(/[ ]{1,},/,",",$0); gsub(/,[ ]{1,}/,",",$0); print}' |tr -s " \+" "\n" |grep -i "^maxDUP=" |perl -pe 's/maxDUP=//i')



####### Parse tool_info.txt file for tool path

TOOL_INFO=${DEDUP_TOOL}"/"tool_info.txt

SAMTOOLS_PATH=$(cat $TOOL_INFO |grep "^SAMTOOLS_PATH" |cut -d '=' -f2)
BEDTOOLS_PATH=$(cat $TOOL_INFO |grep "^BEDTOOLS_PATH" |cut -d '=' -f2)
PICARD_PATH=$(cat $TOOL_INFO |grep "^PICARD_PATH" |cut -d '=' -f2)
JAVA_PATH=$(cat $TOOL_INFO |grep "^JAVA_PATH" |cut -d '=' -f2)
R_PATH=$(cat $TOOL_INFO |grep "^R_PATH" |cut -d '=' -f2)
PERL_PATH=$(cat $TOOL_INFO |grep "^PERL_PATH" |cut -d '=' -f2)



####### remove existing files

if [[ -s ${WORK_DIR}"/"${PEAK}.log.txt ]]
then

rm ${WORK_DIR}"/"${PEAK}.log.txt

fi


if [[ -s ${WORK_DIR}"/"${PEAK}.summary.txt ]]
then

rm ${WORK_DIR}"/"${PEAK}.summary.txt

fi



####### split BAM

echo -e "\n#### Start split BAM ${BAM} into alignments for non-duplicates and those for duplicates, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


${SAMTOOLS}samtools view -h -q ${MAPQ} ${WORK_DIR}"/"${BAM} | \
awk -v f1=${WORK_DIR}"/"${BAM%.bam}.dup.sam -v f2=${WORK_DIR}"/"${BAM%.bam}.nodup.sam 'BEGIN {FS="\t"; OFS="\t"}; \
{if (($1 ~ /^@/) || ($2 ==1024) || ($2==1040)) print $0 >f1; \
if (($1 ~ /^@/) || ($2 ==0) || ($2==16)) print $0 >f2}'


${SAMTOOLS}samtools view -Sbh ${WORK_DIR}"/"${BAM%.bam}.dup.sam > \
${WORK_DIR}"/"${BAM%.bam}.dup.bam

${SAMTOOLS}samtools view -Sbh ${WORK_DIR}"/"${BAM%.bam}.nodup.sam > \
${WORK_DIR}"/"${BAM%.bam}.nodup.bam

rm ${WORK_DIR}"/"${BAM%.bam}.dup.sam
rm ${WORK_DIR}"/"${BAM%.bam}.nodup.sam

echo -e "\nDONE with split BAM ${BAM} into duplicates (.dup.bam) and non-duplicates (.nodup.bam), $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


echo "${BAM}" >> ${WORK_DIR}"/"${PEAK}.summary.txt

${SAMTOOLS}samtools view -c ${WORK_DIR}"/"${BAM%.bam}.nodup.bam | \
awk 'BEGIN {FS=OFS="\t"} {print "Total non-duplicates: "$1}' >> \
${WORK_DIR}"/"${PEAK}.summary.txt

${SAMTOOLS}samtools view -c ${WORK_DIR}"/"${BAM%.bam}.dup.bam | \
awk 'BEGIN {FS=OFS="\t"} {print "Total duplicates: "$1}' >> \
${WORK_DIR}"/"${PEAK}.summary.txt




####### pull out duplicates mapped to peaks

echo -e "\n#### start to extract duplicates from ${BAM%.bam}.dup.bam for peaks ${PEAK}, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


cat ${WORK_DIR}"/"${PEAK} | \
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3}' | \
sort -k1,1 -k2,2n -k3,3n > \
${WORK_DIR}"/"${PEAK}.sorted.bed


cat ${WORK_DIR}"/"${PEAK}.sorted.bed | \
${BEDTOOLS_PATH}"/"intersectBed \
-abam ${WORK_DIR}"/"${BAM%.bam}.dup.bam -b stdin > \
${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam

rm ${WORK_DIR}"/"${BAM%.bam}.dup.bam

echo -e "\nsave the extracted duplicates to ${BAM%.bam}.dup.peak.bam, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

${SAMTOOLS}samtools view -c ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam | \
awk 'BEGIN {FS=OFS="\t"} {print "Total duplicates in peaks: "$1}' >> \
${WORK_DIR}"/"${PEAK}.summary.txt



####### set the the maximum number of duplicates per position at ${maxDUP}
####### if a position has > ${maxDUP} duplicates, the extra duplicates will be removed

echo -e "\n#### start to extract a maximum of ${maxDUP} duplicates per position/strand from ${BAM%.bam}.dup.peak.bam, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

${SAMTOOLS}samtools view ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam | \
awk 'BEGIN {FS=OFS="\t"} {print $3"_"$4"_"$2}' | \
uniq -c | \
awk '{sub(/^[ \t]+/, ""); print}' | \
tr -s " " "\t" | \
awk ' BEGIN {FS="\t"; OFS="\t"} {for (i=1; i<=$1; i=i+1) print i}' > \
${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.dup_freq.txt


${SAMTOOLS}samtools view -H ${WORK_DIR}"/"${BAM} > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.sam


${SAMTOOLS}samtools view ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam | \
paste ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.dup_freq.txt - | \
awk -v f=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.sam -v m=${maxDUP} 'BEGIN {FS=OFS="\t"} {if ($1 <=m) print substr($0, index($0, $2)) >>f}'


${SAMTOOLS}samtools view -Sbh ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.sam > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam 


rm ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.dup_freq.txt
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.sam


echo -e "\nsave the extracted duplicates to ${BAM%.bam}.dup-max${maxDUP}.peak.bam, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


${SAMTOOLS}samtools view -c ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam | \
awk -v m=${maxDUP} 'BEGIN {FS=OFS="\t"} {print "Duplicates in peaks, after setting a max "m" duplicates per position/strand: "$1}' >> \
${WORK_DIR}"/"${PEAK}.summary.txt



####### count number of non-duplicates, duplicates (all), and duplicates by allowing a max of ${maxDUP} duplicates for each peak

echo -e "\n#### start to count #non-duplicates, #all duplicates, and #duplicates (use max ${maxDUP} duplicates per position/strand) for peaks ${PEAK}, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


${SAMTOOLS}samtools view ${WORK_DIR}"/"${BAM%.bam}.nodup.bam | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 ==16) print $3,$4,($4-1+ length($10)),".","1","-"; \
else if ($2 ==0) print $3,$4,($4-1+ length($10)),".","1","+"}' | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 <=0) print $1,"1",$3,$4,$5,$6; else print $0}' > \
${WORK_DIR}"/"${BAM%.bam}.nodup.bam.bed


${SAMTOOLS}samtools view ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 ==1040) print $3,$4,($4-1+ length($10)),".","1","-"; \
else if ($2 ==1024) print $3,$4,($4-1+ length($10)),".","1","+"}' | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 <=0) print $1,"1",$3,$4,$5,$6; else print $0}' > \
${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.bed


${SAMTOOLS}samtools view ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 ==1040) print $3,$4,($4-1+ length($10)),$1,"1","-"; \
else if ($2 ==1024) print $3,$4,($4-1+ length($10)),$1,"1","+"}' | \
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 <=0) print $1,"1",$3,$4,$5,$6; else print $0}' > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.bed



if [[ -s ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt ]]
then

echo -e "\nwarning, ${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt already existed, will delete it, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

rm ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt


elif [[ ! -s ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt ]]
then

echo -e "\nwill save results to ${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

fi




for CHR in $(cat ${WORK_DIR}"/"${PEAK}.sorted.bed |cut -f 1|uniq |tr -s "\n" " ")

do

echo "start $CHR" >> ${WORK_DIR}"/"${PEAK}.log.txt

grep "\b${CHR}\b" ${WORK_DIR}"/"${BAM%.bam}.nodup.bam.bed > \
${WORK_DIR}"/"${BAM%.bam}.nodup.bam.$CHR.bed

grep "\b${CHR}\b" ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.bed > \
${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.$CHR.bed

grep "\b${CHR}\b" ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.bed > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed

cat ${WORK_DIR}"/"${PEAK}.sorted.bed | \
grep "\b${CHR}\b" | \
${BEDTOOLS_PATH}"/"intersectBed -a stdin \
-b ${WORK_DIR}"/"${BAM%.bam}.nodup.bam.$CHR.bed -bed -wa -c | \
${BEDTOOLS_PATH}"/"intersectBed -a stdin \
-b ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.$CHR.bed -bed -wa -c | \
${BEDTOOLS_PATH}"/"intersectBed -a stdin \
-b ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed -bed -wa -c | \
awk 'BEGIN {FS=OFS="\t"} {printf("%s\t%d\t%d\t%.3f\t%d\t%d\t%d\n",$1,$2,$3,($3-$2+1)/1000,$4,$5,$6)}' >> \
${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt

rm ${WORK_DIR}"/"${BAM%.bam}.nodup.bam.$CHR.bed
rm ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.$CHR.bed
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed

done



echo -e "\nDONE, save the #non-duplicates and #duplicates per peak to ${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


rm ${WORK_DIR}"/"${PEAK}.sorted.bed
rm ${WORK_DIR}"/"${BAM%.bam}.nodup.bam.bed
rm ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam.bed
rm ${WORK_DIR}"/"${BAM%.bam}.dup.peak.bam



####### predict the number of duplicates as signal

chmod 755 ${DEDUP_TOOL}"/"lowess.fitting.R

echo -e "\n#### start to predict #duplicates as true signal in peaks using loess function in R, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

cat ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt | \
awk 'BEGIN {FS=OFS="\t"} {printf("%s\t%.2f\t%.2f\n",$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7,(($7/$4)+0.01),(($5/$4)+0.01))}' > \
${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.txt

${R_PATH}"/"R --slave --no-save --no-restore --no-environ --silent --args INPUT=${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.txt OUT=${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.fitted.txt SP=0.2 < ${DEDUP_TOOL}"/"lowess.fitting.R


echo -e "\nDONE, save the prediction to ${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.fitted.txt, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt



####### convert fitted value [duplicates per kb] back to predicted number of duplicates as signal

echo -e "\nreformat the predicted #duplicate-perKb into predicted #duplicate, save to ${PEAK}.duplicate_fitted.txt, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


cat ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.fitted.txt | \
awk 'BEGIN {FS=OFS="\t"} {if(NR>=2) print $0}' | \
awk 'BEGIN {FS=OFS="\t"} {if(NR>=1) gsub(/_/,"\t",$1); print $0}' | \
awk 'BEGIN {FS=OFS="\t"} {printf("%s\t%d\t%d\t%d\t%d\t%d\t%.0f\n",$1,$2,$3,$5,$6,$7,($4*$11))}' | \
awk 'BEGIN {FS=OFS="\t"} {if ($7 >$6) $7=$6; print }' | \
sort -k1,1 -k2,2n -k3,3n | \
awk -v m=${maxDUP} 'BEGIN {FS=OFS="\t"; print "chr\tstart\tend\tnon_duplicate\tduplicate_all\tduplicate_max"m"\tduplicate_predicted_as_signal\tnoise_duplicates"} {print $0,($5-$7)}' > \
${WORK_DIR}"/"${PEAK}.duplicate_fitted.txt


rm ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.fitted.txt 
rm ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count_perkb.txt
rm ${WORK_DIR}"/"${PEAK}.nodup_dup_dup-max${maxDUP}_count.txt



####### get a list of read ID (they are noise duplicates and need to be excluded) 

## get the peaks that had at least one noise duplicate

echo -e "\n#### start to get the read_ID for the noise duplicates that need to be removed from ${BAM%.bam}.dup-max${maxDUP}.peak.bam, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

echo -e "\nfirst extract a list of peaks that contain at least one noise duplicates (i.e., col 6 > col 7) from ${PEAK}.duplicate_fitted.txt, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

cat ${WORK_DIR}"/"${PEAK}.duplicate_fitted.txt | \
awk 'BEGIN {FS=OFS="\t"} {if (NR >=2) print $1,$2,$3,($6-$7)}' | \
awk 'BEGIN {FS=OFS="\t"} {if ($4 >=1) print }' > \
${WORK_DIR}"/"${PEAK}.duplicate_exclude.txt


## intersect these peaks with duplicate alignments

echo -e "\npull out duplicates that overlap the provided peak list, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


for CHR in $(cat ${WORK_DIR}"/"${PEAK}.duplicate_exclude.txt |cut -f 1|uniq |tr -s "\n" " ")

do
echo "$CHR"

grep "\b${CHR}\b" ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.bed | \
${PERL_PATH}perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed

cat ${WORK_DIR}"/"${PEAK}.duplicate_exclude.txt | \
grep "\b${CHR}\b" | \
${BEDTOOLS_PATH}"/"intersectBed -a stdin \
-b ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed -wa -wb -bed >> \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID.bed

rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.$CHR.bed

done


## get a list of read ID for duplicates that are predicted as noise and need to be excluded

awk 'BEGIN {FS=OFS="\t"} {print $1"_"$2"_"$3}' \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID.bed | \
uniq -c | \
awk '{sub(/^[ \t]+/, ""); print}' | \
tr -s " " "\t" | \
awk ' BEGIN {FS="\t"; OFS="\t"} {for (i=1; i<=$1; i=i+1) print i}' | \
paste ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID.bed - | \
awk 'BEGIN {FS=OFS="\t"} {if ($11 <=$4) print $8}' > \
${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID2.bed

echo -e "\nDONE, the read_ID for the noise duplicates is saved to ${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID2.bed, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


rm ${WORK_DIR}"/"${PEAK}.duplicate_exclude.txt
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam.bed
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID.bed




####### exclude noise duplicates from BAM based on read ID

echo -e "\n#### start to excluded noise duplicates from ${BAM%.bam}.dup-max${maxDUP}.peak.bam, based on the list of read_ID in ${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID2.bed, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt

${JAVA_PATH}"/"java \
-jar ${PICARD_PATH}"/"${PICARD}FilterSamReads.jar \
I=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam \
O=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.bam \
EXCLUDE_READS=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID2.bed


rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.reads
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.reads
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_exclude.ID2.bed

${SAMTOOLS}samtools view -c ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.bam | \
awk 'BEGIN {FS=OFS="\t"} {print "Duplicates predicted as signal in peaks: "$1}' >> \
${WORK_DIR}"/"${PEAK}.summary.txt




####### combine the retained duplicates alignments with all the non-duplicates alignments

echo -e "\n#### combine all non-duplicates with duplicates that are predicted to be signals in peaks, $(date)" >> \
${WORK_DIR}"/"${PEAK}.log.txt

echo -e "all non-duplicates are in ${BAM%.bam}.nodup.bam" >> \
${WORK_DIR}"/"${PEAK}.log.txt
echo -e "duplicates predicted as signals are in ${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.bam" >> \
${WORK_DIR}"/"${PEAK}.log.txt


## first sort the retained duplicates alignments

${JAVA_PATH}"/"java \
-jar ${PICARD_PATH}"/"SortSam.jar \
INPUT=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.bam \
OUTPUT=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.sorted.bam \
MAX_RECORDS_IN_RAM=2000000 \
SO=coordinate \
TMP_DIR=$WORK_DIR \
VALIDATION_STRINGENCY=SILENT


## combine the sorted duplicates alignments with the sorted non-duplicates alignments

${JAVA_PATH}"/"java \
-jar ${PICARD_PATH}"/"MergeSamFiles.jar \
I=${WORK_DIR}"/"${BAM%.bam}.nodup.bam \
I=${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.sorted.bam \
O=${WORK_DIR}"/"${BAM%.bam}.deduplicated.bam \
TMP_DIR=$WORK_DIR \
ASSUME_SORTED=TRUE


rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.bam
rm ${WORK_DIR}"/"${BAM%.bam}.nodup.bam
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.bam
rm ${WORK_DIR}"/"${BAM%.bam}.dup-max${maxDUP}.peak.duplicate_kept.sorted.bam

echo -e "\nDONE with the analysis, the output BAM ${BAM%.bam}.deduplicated.bam already coordinates sorted, $(date)\n" >> \
${WORK_DIR}"/"${PEAK}.log.txt


