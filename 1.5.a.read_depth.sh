#!/bin/bash
############################################# Help  ###################################################
Help()
{
   # Display Help
   echo "A very short script to run the samtool command to get depth coverage"
   echo "See 1.5.b.read_depth.R to analyse this output."
   echo "Written by E. Phelps last update 7/9/22"
   echo "~✧ 𓆝 𓆟 𓆞 ~✧"
   echo
   echo "Syntax: ./1.5.a.read_depth.sh [options]"
   echo 
   echo "Options:"
   echo "h     Print this Help."
   echo "s     sample.tab used previously"
   echo "o     Out Directory with path."
   echo "g     reference genome" 
}
###############################################  Options ################################################
# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Semi colon placement mega important if you want to include your own input value.
while getopts ":hg:s:o:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       s) LIST=$OPTARG;;
       o) OUTDIR=$OPTARG;;
       g) GENOME=$OPTARG;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
 done
 
 if [ -z $LIST ] || [ -z $OUTDIR ] ||[ -z $GENOME ]; then 
   echo "Missing options"
   exit
fi
###########################################  Dependents and version ########################################

DATE=$(date "+DATE: %D" | awk '{print $2}')
echo ">1.5.a.read_depth.sh $DATE" >> ${OUTDIR}/resequencing_program_information.out
echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out

#Samtools depth
VER=$(samtools 2>&1 | awk -F"Version:" '{print $2}'| cut -c 1-4| tr -d " \t\n\r" )
echo -e "Samtools\tdepth\t$VER\tcomputes the read depth at each position or region" >> ${OUTDIR}/resequencing_program_information.out

###############################################  Main Program ################################################

if [ -f "depth.out" ]; then
   rm depth.out
fi

echo "Making depth dir"  >> depth.out

if [ ! -d  "${OUTDIR}/depth" ]; then 
      mkdir ${OUTDIR}/depth/
fi 

#See if file exists and has a size greater than zero
echo "Creating a sample list file. Will write to this as we go." >> depth.out
if [ -s "depth.list" ]; then
      rm depth.list
fi

while read LINE; do
   
   SAMPLE_NAME=$(echo $LINE | awk -F" " '{print $1}')
   
   echo ${SAMPLE_NAME}"_depth.txt" >> depth.list
   
   samtools depth -aa ${OUTDIR}/realign/${SAMPLE_NAME}_realigned.bam --reference ${GENOME} > ${SAMPLE_NAME}'_depth.txt'
   
done < ${LIST}

mv *depth.txt ${OUTDIR}/depth/

echo "Finished calculating depth... Move to Rscript 1.5b. ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> depth.out

