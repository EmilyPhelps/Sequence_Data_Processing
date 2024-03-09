#!/bin/bash
############################################# Help  ###################################################
# running test in screen -r 36622.pts-6.bio-18mit
Help()
{
   # Display Help
   echo "Mapping script for low depth resequencing projects in the Taylor lab 2022."
   echo " Script trims and mapps the sequences "
   echo "Written by E.Phelps, last update 14/08/2023"
   echo "~✧ 𓆝 𓆟 𓆞 ~✧"
   echo
   echo "Syntax: ./1.1.mapping.sh [options]"
   echo 
   echo "Options:"
   echo "h     Print this Help."
   echo "g     Genome, with path!"
   echo "s     List with each line being a separate sample. Line must match"
   echo "      sample name exactly. Must be in a folder with the matching title."
   echo "b     Batch number- a user input."
   echo "l     Library name"
   echo "o     Out Directory with path."
   echo "t     Temp directory"
}
###############################################  Options ################################################

# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Colon placement mega important if you want to include your own input value.
while getopts ":hg:s:b:l:o:t:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       g) GENOME=$OPTARG;;
       s) LIST=$OPTARG;;
       b) BATCH=$OPTARG;;
       l) LIB=$OPTARG;;
       o) OUTDIR=$OPTARG;;
       t) TEMP_DIR=$OPTARG;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
done

if [ -z $GENOME ] || [ -z $LIST ] || [ -z $OUTDIR ] || [ -z $TEMP_DIR ] ; then 
   echo "Missing options"
   exit
fi

###########################################  Dependents #################################################

PATH=/usr/bin/:$PATH
DATE=$(date "+DATE: %D" | awk '{print $2}')
echo ">1.1.mapping.sh $DATE" >> ${OUTDIR}/resequencing_program_information.out
echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out
#JAVA
VER=$(java -version 2>&1 | awk -F '"' 'NR==1 {print $2}')
echo -e "Java\t-\t$VER\t-" >> ${OUTDIR}/resequencing_program_information.out

#BWA indexing 
#tr to remove white space
VER=$(bwa 2>&1 | awk -F"Version:" '{print $2}'| tr -d " \t\n\r" )
echo -e "BWA\tidex\t$VER\tIndex ref genome" >> ${OUTDIR}/resequencing_program_information.out
#BWA mem
echo -e "BWA\tmem\t$VER\tMap sequences to ref genome" >> ${OUTDIR}/resequencing_program_information.out

#PICARD CreateSequenceDictionary
VER=$(java -Xmx8G -jar /usr/bin/picard.jar CreateSequenceDictionary 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}' | tr -d " \t\n\r" )
echo -e "PICARD\tCreateSequenceDictionary\t$VER\tCreates a sequence dictionary for a reference sequence" >> ${OUTDIR}/resequencing_program_information.out

#PICARD FastqToSam
VER=$(java -Xmx8G -jar /usr/bin/picard.jar FastqToSam 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}' | tr -d " \t\n\r" )
echo -e "PICARD\tFastqToSam\t$VER\tConverts a FASTQ file to an unaligned BAM or SAM file" >> ${OUTDIR}/resequencing_program_information.out

#PICARD RevertSam
VER=$(java -Xmx8G -jar /usr/bin/picard.jar RevertSam 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}'| tr -d " \t\n\r" )
echo -e "PICARD\tRevertSam\t$VER\treverts sam to bam file, removing some information" >> ${OUTDIR}/resequencing_program_information.out

#PICARD MarkIlluminaAdapters
VER=$(java -Xmx8G -jar /usr/bin/picard.jar MarkIlluminaAdapters 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}' | tr -d " \t\n\r" )
echo -e "PICARD\tMarkIlluminaAdapters\t$VER\tMarks adapters" >> ${OUTDIR}/resequencing_program_information.out

#PICARD SamToFastq
VER=$(java -Xmx8G -jar /usr/bin/picard.jar SamToFastq 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}' | tr -d " \t\n\r" )
echo -e "PICARD\tSamToFastq\t$VER\tCoverts a Sam file to a Fastq file - so it can be aligned" >> ${OUTDIR}/resequencing_program_information.out

#PICARD MergeBamAlignment
VER=$(java -Xmx8G -jar /usr/bin/picard.jar MergeBamAlignment 2>&1 | awk -F"Version:" '{print $2}' | awk -F" " '{print $1}'| tr -d " \t\n\r" )
echo -e "PICARD\tMergeBamAlignment\t$VER\t Merge alignment data from a SAM or BAM with data in an unmapped BAM file" >> ${OUTDIR}/resequencing_program_information.out
###########################################  Main program #################################################

#1. make dict file for reference genome - only needs to be done once!
if [ -f "mapping.out" ]; then
   rm mapping.out
fi

echo "Checking for genome dictionary." >> mapping.out


if [[ -f "${GENOME}.*" ]] && [ -f ${GENOME}.amb ];then
   echo "Index exists" >> mapping.out
else
   echo "indexing genome"
   bwa index $GENOME
fi

#The original code here created a file with a .fa.dict suffix which was not picked up. 
#Must be just .dict
G_NAME=$(echo $GENOME| awk -F".fa" '{print $1}')
if [ -f ${G_NAME}.dict ]; then 
 if [ -s ${G_NAME}.dict ]; then
   echo "Directory present- moving on" >> mapping.out
 else
   rm ${G_NAME}.dict 
   echo "No genome dictionary found. Making genome dictionary"
   java -jar /usr/bin/picard.jar CreateSequenceDictionary \
       R=${GENOME} \
       O=${G_NAME}.dict
 fi
else
 echo "No genome dictionary found. Making genome dictionary" >> mapping.out
 java -jar /usr/bin/picard.jar CreateSequenceDictionary \
       R=${GENOME} \
       O=${G_NAME}.dict
fi

for SAMPLE in `cat $LIST`; do 
   SAMPLE_NAME=$(echo $SAMPLE | awk -F"/" '{print $(NF)}')
   echo "Starting with ${SAMPLE}" >> mapping.out
  
#2. Collate sample information and create out directories
   #BATCH, LANE, SAMPLE NAME, ETC..
   #check lane code!

   echo "Creating sample_information.tab file" >> mapping.out
   
   if [ -z $BATCH ]; then
      BATCH=NA
   fi
   
   if [ -z $LIB  ]; then
      LIB=NA
   fi
   
   LANE=$(ls ${SAMPLE}/*fq.gz | tail -c 11 | cut -c1-2)
   
   if [ -z sample_information.tab ]; then
     echo -e "SAMPLE_NAME\tBATCH\tLANE\tLIBRARY" >> sample_information.tab
     echo -e ${SAMPLE_NAME}"\t"$BATCH"\t"$LANE"\t"$LIB >> sample_information.tab
    else
     echo -e ${SAMPLE_NAME}"\t"$BATCH"\t"$LANE"\t"$LIB >> sample_information.tab
   fi

   echo -e "Creating output directory- ${OUTDIR}." >> mapping.out

   if [ -d $OUTDIR ]; then
    echo
    if [ "$(ls -A $OUTDIR)" ]; then
       echo -e "$OUTDIR exists and is not empty. Do what you want with this info" >> mapping.out
    else 
       echo -e "Directory exists but is empty" >> mapping.out
     fi
   else
      mkdir ${OUTDIR}
   fi
 
#3. Convert to sam format and remove problematic records

   echo "Converting fastq to a bam files" >> mapping.out
   if [ ! -d  "${OUTDIR}/fastq_bam" ]; then
      mkdir ${OUTDIR}/fastq_bam
   fi
   
   if [ ! -d "${OUTDIR}/fastq_bam/${LIB}" ]; then
      mkdir ${OUTDIR}/fastq_bam/${LIB}
   fi
 

  LANE_ARR=( $(find ${SAMPLE}/*1.fq.gz | sed 's/_1.fq.gz//g'))
  
   for FILE in ${LANE_ARR[@]}; do
      F1=$(echo ${FILE}_1.fq.gz)
      F2=$(echo ${FILE}_2.fq.gz)
      LANE=$(echo ${FILE} | awk '{ print substr( $0, length() - 1) }')
     java -Xmx8G -jar /usr/bin/picard.jar FastqToSam \
         FASTQ=${F1} \
         FASTQ2=${F2} \
         OUTPUT=./${SAMPLE_NAME}_${LANE}_fastqtosam.bam \
         SAMPLE_NAME=${SAMPLE_NAME} \
         LIBRARY_NAME=${LIB} \
         PLATFORM=illumina \
         SEQUENCING_CENTER=Novogene \
         RUN_DATE=2022-06-26T00:00:00-0400 \
         TMP_DIR=$TEMP_DIR

      mv *_fastqtosam.bam ${OUTDIR}/fastq_bam/${LIB}/

      if [ ! -d  "${OUTDIR}/revert_sam" ]; then 
         mkdir ${OUTDIR}/revert_sam
      fi

      if [ ! -d  "${OUTDIR}/revert_sam/${LIB}" ]; then
         mkdir ${OUTDIR}/revert_sam/${LIB}
      fi

      echo "Converting BAM to uBAM and remove problematic records via RevertSam" >> mapping.out

      java -Xmx8G -jar /usr/bin/picard.jar RevertSam \
          I=${OUTDIR}/fastq_bam/${LIB}/${SAMPLE_NAME}_${LANE}_fastqtosam.bam \
          O=${SAMPLE_NAME}_${LANE}_revertsam.bam \
          SANITIZE=true \
          MAX_DISCARD_FRACTION=0.005 \
          ATTRIBUTE_TO_CLEAR=XT \
          ATTRIBUTE_TO_CLEAR=XN \
          ATTRIBUTE_TO_CLEAR=AS \
          ATTRIBUTE_TO_CLEAR=OC \
          ATTRIBUTE_TO_CLEAR=OP \
          SORT_ORDER=queryname \
          RESTORE_ORIGINAL_QUALITIES=true \
          REMOVE_DUPLICATE_INFORMATION=true \
          REMOVE_ALIGNMENT_INFORMATION=true \
          TMP_DIR=$TEMP_DIR

      mv *.bam ${OUTDIR}/revert_sam/${LIB}/

      # 4. Remove adapters
      echo "Removing adapters" >> mapping.out

     if [ ! -d  "${OUTDIR}/mark_adapt" ]; then 
         mkdir ${OUTDIR}/mark_adapt
     fi

     if [ ! -d   "${OUTDIR}/mark_adapt/${LIB}" ]; then
         mkdir ${OUTDIR}/mark_adapt/${LIB}
     fi 

      java -Xmx8G -jar /usr/bin/picard.jar MarkIlluminaAdapters \
         I=${OUTDIR}/revert_sam/${LIB}/${SAMPLE_NAME}_${LANE}_revertsam.bam  \
         O=${SAMPLE_NAME}_${LANE}_mark_illumina_adapters.bam \
         ADAPTERS=null \
         FIVE_PRIME_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
         THREE_PRIME_ADAPTER=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
         M=${SAMPLE_NAME}_mark_illumina_adapters.txt \
         TMP_DIR=$TEMP_DIR

      mv *_mark_illumina_adapters.bam ${OUTDIR}/mark_adapt/${LIB}


     #Converting BAM to FASTQ and removing adapter sequencing using SamToFastq is the first part of the piped command below.

    # 4. Map reads using bwa
    echo "Lets map some reads" >> mapping.out
    set -o pipefail
    #this is a piped command
      java -Xmx8G -jar /usr/bin/picard.jar SamToFastq \
       I= ${OUTDIR}/mark_adapt/${LIB}/${SAMPLE_NAME}_${LANE}_mark_illumina_adapters.bam\
       FASTQ=/dev/stdout \
      CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
       TMP_DIR=$TEMP_DIR | \
       
      bwa mem -M -t 63 -p $GENOME /dev/stdin | \

      java -Xmx16G -jar /usr/bin/picard.jar MergeBamAlignment \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${OUTDIR}/revert_sam/${LIB}/${SAMPLE_NAME}_${LANE}_revertsam.bam \
        OUTPUT=${SAMPLE_NAME}_${LANE}_piped.bam \
        R=$GENOME CREATE_INDEX=true ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
        TMP_DIR=$TEMP_DIR

     if [ ! -d  "${OUTDIR}/piped_bam" ]; then 
         mkdir ${OUTDIR}/piped_bam
     fi 

     if [ ! -d  "${OUTDIR}/piped_bam/${LIB}" ]; then 
         mkdir ${OUTDIR}/piped_bam/${LIB}
     fi 

     mv *_piped.bam ${OUTDIR}/piped_bam/${LIB}/

     echo "Completed mapping sample ${SAMPLE_NAME}_${LANE}" >> mapping.out
 done 
done

echo "~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ Mapping complete ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> mapping.out
