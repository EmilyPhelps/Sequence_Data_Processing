#!/bin/bash
############################################# Help  ###################################################

Help()
{
   # Display Help
   echo "This script identifies multiple bam files and merges them."
   echo "These files occur when aliquots of a sample have been run on multiple lanes"
   echo "Files must be merged before duplicate removal"
   echo "~ 𓆜 𓆜<3~"
   echo
   echo "Syntax: 1.1.a.merge.sh "
   echo 
   echo options
   echo "h     Print this Help."
   echo "b     Path to folder containing bams that need to be merged"
   
   echo
   echo ""
}
###############################################  Options ################################################

# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Colon placement mega important if you want to include your own input value.
while getopts ":hb:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       b) DIR=$OPTARG
       ;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
done

if [ -z $DIR ]; then 
   echo "Missing bam directory"
   exit
fi
#######################################  Dependencies and Version. ##########################################

DATE=$(date "+DATE: %D" | awk '{print $2}')

echo ">1.1.a.merge.sh $DATE" >> ${OUTDIR}/resequencing_program_information.out

echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out

#Samtools merge
VER=$(samtools 2>&1 | awk -F"Version:" '{print $2}'| cut -c 1-4 | tr -d " \t\n\r" )
echo -e "Samtools\tmerge\t$VER\tmerge multi-lane bams" >> ${OUTDIR}/resequencing_program_information.out

###########################################  Main program. #################################################
rm merge.out

echo "Begining merging" >> merge.out

SAMPLE_ARR=( $(ls ${DIR}*bam | awk -F"/" '{print $(NF)}' | awk -F"_" '{print $1}' |sort | uniq) )

for SAMPLE in ${SAMPLE_ARR[@]} ; do
   FIND=( $(find ${DIR} -name ${SAMPLE}_*bam | awk -F"/" '{print $(NF)}') )
   if [ "${#FIND[@]}" -gt 1 ]; then 
     #Length of the array minus 1
     NUM=$((${#FIND[@]}-1))
     for i in $(seq 0 $NUM); do
      echo ${DIR}/${FIND[${i}]} >> merge.list.tmp
     done
     echo "merging "${SAMPLE} >> merge.out
     samtools merge -b merge.list.tmp ${SAMPLE}_piped.bam
     echo "complete "${SAMPLE} >> merge.out
     rm merge.list.tmp
   fi
done   

echo "~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ Merging complete ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> merge.out
      

