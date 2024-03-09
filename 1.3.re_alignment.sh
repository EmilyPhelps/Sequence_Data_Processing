#!/bin/bash
############################################# Help  ###################################################

Help()
{
   # Display Help
   echo "Written by E.Phelps, last update 24/08/2022"
   echo "This script requires an older version of GATK "
   echo "(GATK3) and also Java (Java 8). If you get the"
   echo "error relating to not being able to find the Indel"
   echo "modules this is likely the cause. See Dependents part of "
   echo "script"
   echo "~✧ 𓆝 𓆟 𓆞 ~✧"
   echo
   echo "Syntax: ./1.3.re_alignment.sh [options]"
   echo 
   echo "Options:"
   echo "h     Print this Help."
   echo "s     This list is different! It is the sample_information.tab from mapping."
   echo "      This is so you can submit all the samples in one command now"
   echo "o     Out Directory"
   echo "t     Temp directory"
   echo "g      obviously genome"
}
###############################################  Options ################################################
# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Colon placement mega important if you want to include your own input value.
while getopts ":hg:s:o:t:g:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       s) LIST=$OPTARG;;
       o) OUTDIR=$OPTARG;;
       t) TEMP_DIR=$OPTARG;;
       g) GENOME=$OPTARG;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
done

if [ -z $LIST ] || [ -z $OUTDIR ]|| [ -z $TEMP_DIR ]|| [ -z $GENOME ]; then 
   echo "Missing options"
   exit
fi

#########################################  Dependents and versions ######################################
#Replace these with the paths to your programs. 

PATH=/usr/bin/:$PATH
GATK3=/local/fatbob/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar 
JAVA8=~/java8/jre1.8.0_321/bin/java

DATE=$(date "+DATE: %D" | awk '{print $2}')
echo ">1.3.re_alignment.sh $DATE" >> ${OUTDIR}/resequencing_program_information.out
echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out

#JAVA
VER=$(~/java8/jre1.8.0_321/bin/java -version 2>&1 | awk -F '"' 'NR==1 {print $2}')
echo -e "Java\t-\t$VER\t-" >> ${OUTDIR}/resequencing_program_information.out

#Samtools faidx
VER=$(samtools 2>&1 | awk -F"Version:" '{print $2}'| cut -c 1-4 | tr -d " \t\n\r" )
echo -e "Samtools\tfaidx\t$VER\tIndex reference" >> ${OUTDIR}/resequencing_program_information.out

#GATK3 RealignerTargetCreator
VER=$(${JAVA8} -jar ${GATK3} -v 2>&1 | awk -F"version" '{print $2}' | sed 's/)://g'| tr -d " \t\n\r" )
echo -e "GATK3\tRealignerTargetCreator\t$VER\tidentify indel regions" >>  ${OUTDIR}/resequencing_program_information.out

#GATK3 IndelRealigner
echo -e "GATK3\tIndelRealigner\t$VER\tRealign around indels" >>  ${OUTDIR}/resequencing_program_information.out

######################################### Main Program ###################################################
if [ -f "realign.out" ]; then
   rm realign.out
fi

#1. Index reference 
echo "Indexing reference" >> realign.out

if [ ! -f ${GENOME}.fai ]; then
   samtools faidx ${GENOME}
fi

echo "Creating a sample list" >> realign.out

while read LINE; do
   SAMPLE_NAME=$(echo $LINE | awk -F" " '{print $1}')
   LIB=$(echo $LINE | awk -F" " '{print $4}')
   
   #2. Index bam
   samtools index ${OUTDIR}/clipped/${LIB}/${SAMPLE_NAME}_clipped.bam
 
   #3. Create a sample_list to feed into GATK
   echo ${OUTDIR}/clipped/${LIB}/${SAMPLE_NAME}_clipped.bam >> sample.list   
  
done < ${LIST}

echo "Creating a list of indels" >> realign.out

#4.Create a list of indels

${JAVA8} -jar ${GATK3} -T RealignerTargetCreator \
	-R ${GENOME} \
	-I sample.list \
	-o ${OUTDIR}/indel_realigner.intervals

echo "Completed indel list" >> realign.out	

echo "Re-aligning indels" >> realign.out

#5.Re-align

${JAVA8} -jar ${GATK3}  -T IndelRealigner \
	-R ${GENOME} \
	-I sample.list \
	-targetIntervals ${OUTDIR}/indel_realigner.intervals\
	--consensusDeterminationModel USE_READS\
	--nWayOut _realigned.bam

echo "Completed re-alignment" >> realign.out	
echo "Just doing a little house keeping" >> realign.out

#6.Rename and organise
if [ ! -d  "${OUTDIR}/realign" ]; then 
      mkdir ${OUTDIR}/realign/
fi 

mv *realigned.ba* ${OUTDIR}/realign/

SAMPLE_ARR=( $(find ${OUTDIR}/realign/ -name '*ba*' | awk -F"/" '{print $(NF)}') )
  
  for BA in ${SAMPLE_ARR[@]}; do
  	NEW_NAME=$(echo $BA | sed 's/_clipped//g' )
	mv $BA $NEW_NAME
  done	


rm sample.list
echo "~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ Re-alignment Complete ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> realign.out

