#!/bin/bash
############################################# Help  ###################################################
Help()
{
   # Display Help
   echo "Unfortunately you cant run lacer.pl from outside the docker..."
   echo "This script writes a mega long script that runs lacer separately for each file"
   echo "It a bit convoluted because it was written to be executed from outside docker"
   echo "Once this script is complete enter docker and run ./home/recal_subscript.sh"
   echo "~✧ 𓆝 𓆟 𓆞 ~✧"
   echo
   echo "Syntax: ./1.4.recalibration.sh [options]"
   echo 
   echo "Options:"
   echo "h     Print this Help."
   echo "o     Out Directory- sensitive so dont have a slash at the end"
   echo "g     Genome"
}
###############################################  Options ################################################
# Add options in by adding another letter eg. :x, and then include a letter) and some commands.
# Colon placement mega important if you want to include your own input value.
while getopts ":ho:g:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
       o) OUTDIR=$OPTARG
       ;;
       g) GENOME=$OPTARG
       ;;
      \?) # incorrect option
         echo "Error: Invalid option щ(ಥДಥщ)"
         exit;;
   esac
done

if [ -z ${OUTDIR} ]|| [ -z ${GENOME} ]; then
   echo "Missing options"
   exit 1
fi

###########################################  Dependents ##################################################
#This script requires DOCKER!

DATE=$(date "+DATE: %D" | awk '{print $2}')
echo ">1.4.recalibration $DATE" >> ${OUTDIR}/resequencing_program_information.out
echo -e "Software\tProgram\tVersion\tRole" >> ${OUTDIR}/resequencing_program_information.out

#Docker
VER=$(docker -v | awk -F" " '{print $3}' | sed 's/,//g')
echo -e "docker\t-\t${VER}\t-" >> ${OUTDIR}/resequencing_program_information.out

#Lacer- No inprogram version info
VER="0.426"
echo -e "Lacer\t-\t${VER}\tCalculate Recalibration scores and create table" >> ${OUTDIR}/resequencing_program_information.out

#Lacepr-No inprogram version info
VER="0.2"
echo -e "Lacepr\t-\t${VER}\tIntegrate new scores with bam files" >> ${OUTDIR}/resequencing_program_information.out
###############################################  Main Program ############################################
#Include dir and out file organisation
if [ ! -d "${OUTDIR}/realign" ]; then
   echo "Have you run re-align? No directory present"
   exit 1
fi

if [ -f "recalibration.out" ]; then
   rm recalibration.out
fi

#1. Set up docker container. Snooze.

echo "Setting up docker container" >> recalibration.out

if [ "$(docker ps -a | grep ${OUTDIR})" ]; then
  if [ "$(docker ps -aq -f status=exited -f name=${OUTDIR})" ]; then
   docker rm ${OUTDIR}
  else
   docker stop ${OUTDIR}
   docker rm ${OUTDIR}
  fi
fi

docker run -v ${PWD}/${OUTDIR}/realign:/home/realign --name ${OUTDIR} -d -it andreaswilm/lacer > cont_id

CONT_ID=$(cat cont_id)

docker cp ${GENOME} ${CONT_ID}:/home

if [ ! -f ${GENOME}.fai ]; then
 echo "No genome index, re-indexing" >> recalibration.out
 samtools faidx ${GENOME}
fi

docker cp ${GENOME}.fai ${CONT_ID}:/home

#2. Create a sub script containing a for loop to upload and run in docker
echo "writing sub-script for to run in docker" >> recalibration.out

SAMPLE_ARR=( $(find ${OUTDIR}/realign/ -name '*bam' | awk -F"/" '{print $(NF)}') )

if [ -f recal_subscript.sh ]; then
   rm recal_subscript.sh
fi

echo '#!/bin/bash' >> recal_subscript.sh
echo -e "#1. Calibraton table" >> recal_subscript.sh

for BAM in ${SAMPLE_ARR[@]}; do
   SAMPLE_NAME=${BAM::-4}
   echo "#creating ${BAM} calibration table" >> recal_subscript.sh
   echo -e "lacer.pl -bam realign/${BAM} -reference ${GENOME} -outfile ${SAMPLE_NAME}.txt" >> recal_subscript.sh
   echo -e "if [ -s ${SAMPLE_NAME}.txt ]; then\n echo \"lacer.pl failure\"\n exit 1\nfi\n" >> recal_subscript.sh
done

echo -e "\n#2. Recalibrate\n" >> recal_subscript.sh

for BAM in ${SAMPLE_ARR[@]}; do
   SAMPLE_NAME=${BAM::-4}
   echo "#recalibrating ${BAM}" >> recal_subscript.sh
   echo -e "lacepr --bam realign/${BAM} --recal ${SAMPLE_NAME}.txt --out ${SAMPLE_NAME}_recal.bam\n" >> recal_subscript.sh
done

echo -e "\n#3.House keeping" >> recal_subscript.sh
echo -e "mkdir recalibration\nmv *recal.bam recalibration/" >> recal_subscript.sh

chmod 755 recal_subscript.sh 

docker cp recal_subscript.sh ${CONT_ID}:/home

docker exec -d ${CONT_ID} chmod /home/recal_subscript.sh 755

#3. Executing script and recovering files
#echo "Starting to create recalibration tables and recalibrating samples"  >> recalibration.out
# docker exec -d ${OUTDIR} ./home/recal_subscript.sh

#echo "Finished recalibration- Housekeeping time"  >> recalibration.out

#if [ ! -d  "${OUTDIR}/recalibration" ]; then 
#      mkdir ${OUTDIR}/recalibration
#fi 
 
#docker cp ${CONT_ID}:/home/recalibration ${OUTDIR}/recalibration/

#rm cont_id
#rm recal_subscript.sh
echo "~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ Recalibration prep complete ~✧ 𓆝 𓆟 𓆞 𓆝 𓆟 𓆞 ~✧ " >> recalibration.out
