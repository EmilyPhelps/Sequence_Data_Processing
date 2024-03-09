#!/bin/bash
############################################# Help  ###################################################
Help()
{
   # Display Help
   echo "This script filters the poorly aligned sequences and removes duplicates"
   echo "This filters out alignments with a MAPQ score lower than 20. This is the"
   echo "equivalent of correct mapping probability 0.99. If BamUtil stage fails make sure the path"
   echo "for bamUtil is correct. Always check the size of the piped bam. Some may fail if there are"
   echo "multiple lane runs e.g. CC7 and CC8 in the Fulleri Proj"
   echo "Written by E.Phelps, last update 05/09/2022"
   echo "~✧ 𓆝 𓆟 𓆞 ~✧"
   echo
   echo "Syntax: ./1.2.duplicates.sh [options]"
   echo 
   echo "Options:"
   echo "h     Print this Help."
   echo "s     List with each line being a separate sample. Line must match"
   echo "      sample name exactly. Must be in a folder with the matching title."
   echo "l     Library name"
   echo "o     Out Directory with path."
   echo "t     Temp directory"
}
###############################################  Options ################################################
# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Colon placement mega important if you want to include your own input value.
while getopts ":hg:l:s:o:t:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       s) LIST=$OPTARG;;
       o) OUTDIR=$OPTARG;;
       l) LIB=$OPTARG;;
       t) TEMP_DIR=$OPTARG;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
done

if [ -z $LIST ] || [ -z $OUTDIR ] || [ -z $TEMP_DIR ]|| [ -z $LIB ]; then 
   echo "Missing options"
   exit
fi

###########################################  Dependents and version ########################################

PATH=/usr/bin/:$PATH
BAM_UTIL=/local/chh20csu/bamUtil/bin/bam

DATE=$(date "+DATE: %D" | awk '{print $2}')
echo ">1.2.duplicates.sh $DATE" >> ${OUTDIR}/resequencing_program_information.out
echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out

#Samtools view
VER=$(samtools 2>&1 | awk -F"Version:" '{print $2}'| cut -c 1-4 | tr -d " \t\n\r" )
echo -e "Samtools\tview\t$VER\tfiltering bams" >> ${OUTDIR}/resequencing_program_information.out

#Java
VER=$(java -version 2>&1 | awk -F '"' 'NR==1 {print $2}'| tr -d " \t\n\r" )
echo -e "Java\t-\t$VER\t-" >> ${OUTDIR}/resequencing_program_information.out

#Picard MarkDuplicates
VER=$(java -Xmx8G -jar /usr/bin/picard.jar MarkDuplicates 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}' |tr -d " \t\n\r" )
echo -e "PICARD\tMarkDuplicates\t$VER\tIdentifies and Marks Duplicates" >> ${OUTDIR}/resequencing_program_information.out

#BamUtil ClipOverlap
VER=$($BAM_UTIL 2>&1 | awk -F"Version:" '{print $2}'| awk -F";" '{print $1}'| tr -d " \t\n\r" )
echo -e "BamUtil\tClipOverlap\t$VER\tClips overlapping read pars in BAM files" >> ${OUTDIR}/resequencing_program_information.out

###########################################  Main program #################################################

if [ -f "duplicates.out" ]; then
   rm duplicates.out
fi

for SAMPLE in `cat $LIST`; do
 #1. Filter the bam files
 # -h, -u change the way the sam file is viewed. -q is remove those with map score lower than 20. -S ignore compatibility
 # with older samtools. 
 
 echo "Starting with ${SAMPLE}" >> duplicates.out
 
 echo "Filtering..." >> duplicates.out
 
 SAMPLE_NAME=$(echo $SAMPLE | awk -F"/" '{print $(NF)}')
 BAM=${OUTDIR}/piped_bam/${LIB}/${SAMPLE_NAME}*_piped.bam
 
 if [ ! -d  "${OUTDIR}/piped_bam" ] && [ -z ${BAM} ]; then 
      echo "piped file does not exist- have you run 1.1.mapping?"
      exit
 fi 
 
 if [ ! -d  "${OUTDIR}/filt_bam" ]; then 
      mkdir ${OUTDIR}/filt_bam/
 fi 
 
 if [ ! -d  "${OUTDIR}/filt_bam/${LIB}" ]; then 
      mkdir ${OUTDIR}/filt_bam/${LIB}
 fi 
 
 time samtools view -h -q 20 ${BAM} | samtools view -buS | samtools sort -o ${SAMPLE_NAME}_minq20_sorted.bam
 
 mv *_minq20_sorted.bam ${OUTDIR}/filt_bam/${LIB}/
 
 echo "Filtering...completed! " >> duplicates.out
 
 #2. Mark duplicates
 
 echo "Removing duplicates" >> duplicates.out
 
 if [ ! -d  "${OUTDIR}/remove_dups" ]; then 
      mkdir ${OUTDIR}/remove_dups/
 fi 
 
 if [ ! -d  "${OUTDIR}/remove_dups/${LIB}" ]; then 
      mkdir ${OUTDIR}/remove_dups/${LIB}
 fi 
 
 #Removes both optical and sequencing duplicates
 
 java -Xmx8G -jar /usr/bin/picard.jar MarkDuplicates \
   I=${OUTDIR}/filt_bam/${LIB}/${SAMPLE_NAME}_minq20_sorted.bam \
   O=${SAMPLE_NAME}_remove_dups.bam \
   M=${SAMPLE_NAME}_remove_dups_stat.txt \
   VALIDATION_STRINGENCY=SILENT \
   REMOVE_DUPLICATES=true

   mv *remove_dups* ${OUTDIR}/remove_dups/${LIB}/
  
   echo "Duplicates removed" >> duplicates.out
   echo "Clipping reads"
   
   if [ ! -d  "${OUTDIR}/clipped" ]; then 
      mkdir ${OUTDIR}/clipped/
   fi 
 
   if [ ! -d  "${OUTDIR}/clipped/${LIB}" ]; then 
      mkdir ${OUTDIR}/clipped/${LIB}
   fi 
   
   $BAM_UTIL clipOverlap \
      --in ${OUTDIR}/remove_dups/${LIB}/${SAMPLE_NAME}_remove_dups.bam \
      --out ${OUTDIR}/clipped/${LIB}/${SAMPLE_NAME}_clipped.bam\
      --stats
   
   mv *clipped* ${OUTDIR}/clipped/${LIB}/
   echo "Clipped  <3"
   
   echo "${SAMPLE_NAME} processed" >> duplicates.out
 done
 
echo "~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ Filtering and Duplication removal complete ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> duplicates.out
