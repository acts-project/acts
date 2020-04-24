#!/bin/bash
#
# This script runs the KalmanFitter timing test with different smoothing method at different |eta| or pT bins
# using particle gun particles
#
helpFunction()
{
   echo ""
   echo "Usage: $0 -v <binVal> -d <detector> -b <bFieldMap> -t <numTracksPerEvent> -n <numEvents>"
   echo -e "\t-v The binning variable, either 'eta' or 'pt'. Required"
   echo -e "\t-d The detector type, either 'Generic' or 'DD4hep'. Required" 
   echo -e "\t-b The '.txt' or '.root' file for B Field map. Optional. In default using constant BField: (0, 0, 2)"
   echo -e "\t-t The number of tracks per event. Optional. In default: 100"
   echo -e "\t-n The number of events. Optional. In default: 1"
   exit 1 # Exit script after printing help
}

#Default number of tracks/event and event
numTracksPerEvent=100
numEvents=1

# Get parameters
while getopts "v:d:b:t:n" opt
do
   case "$opt" in
      v ) binVal="$OPTARG" ;;
      d ) detector="$OPTARG" ;;
      b ) bFieldMap="$OPTARG" ;;
      t ) numTracksPerEvent="$OPTARG" ;;
      n ) numEvents="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case unknown parameter
   esac
done

#Check input for binning variable
if [ -z "${binVal}" ]; then
   echo "binVal is empty";
   helpFunction
elif [[ ${binVal} -ne "eta" ]] && [[ ${binVal} -ne "pt" ]]; then
   echo "binVal is not valid. It should be either 'eta' or 'pt'!"
   exit 1
fi

# Check input for detector 
if [ -z "${detector}" ]; then
   echo "detector is empty";
   helpFunction
fi

# Check input for B field
bField='--bf-value 0 0 2'
if [ -z "${bFieldMap}" ]; then
   echo "bFieldMap is empty. Will use constant B Field (0, 0, 2) then."
elif [ -f "$bFieldMap" ] ; then
   echo "Input bField map: $bFieldMap "
   bField='--bf-map ${bFieldMap}'
else
   echo "Input file for bField map file does not exist!"
   exit 1
fi

#Print the parameters 
echo "Binning Variable: ${binVal}"
echo "Test detector: ${detector}"
echo "B Field: ${bField}"
echo "Number of tracks per event: ${numTracksPerEvent}"
echo "Number of events: ${numEvents}"

time_stamp=`date +%F%T`
exe_directory=$PWD
run_directory=KalmanFitter_timing_${time_stamp}
mkdir ${run_directory}
cd ${run_directory}

output_file='output.log'
echo "*****KalmanFitter timing test vs. ${binVal}*****" > ${output_file}
echo "* Test Detector: ${detector}" >> ${output_file}
echo "* B Field: ${bField}" >> ${output_file}
echo "* Events: ${numEvents}" >> ${output_file}
echo "* Tracks/event: ${numTracksPerEvent}" >> ${output_file} 
echo "****************************************" >> ${output_file}
echo "*"
echo "* job | ${bFieldMap} " >> ${output_file}

jobID=0
if [[ ${binVal} -eq "eta" ]];then
   # Loop over the eta bins 
   for eta in 0 0.5 1 1.5 2 2.5 ; do
      # Run simulation
      simulation="${exe_directory}/ActsSimFatras${detector} ${bField} -n ${numEvents} --evg-input-type gun --pg-nparticles ${numTracksPerEvent} --pg-pt-range 0.1 100 --pg-eta-range ${eta} ${eta}  --fatras-pmin-gev 0.1 --output-csv=1 --output-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}"  
      echo ${simulation}
      eval ${simulation}
      # Run reconstruction 
      reconstruction="$exe_directory/ActsRecTruthTracks ${bField} --input-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta} --output-dir reco_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}"
      eval ${reconstruction}
      
      echo "${jobID}, ${eta}" >> ${output_file}
      # Archive with Job ID
      mv reco_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}/timing.tsv timing_${jobID}.tsv
      # JobID
      let "jobID++"
   done
else 
   # Loop over the pt bins 
   for pt in 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 50.0 100.0 ; do
      # Run simulation
      simulation="${exe_directory}/ActsSimFatras${detector} ${bField} -n ${numEvents} --evg-input-type gun --pg-nparticles ${numTracksPerEvent} --pg-pt-range ${pt} ${pt} --pg-eta-range -2.5 2.5  --fatras-pmin-gev 0.1 --output-csv=1 --output-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_pt${pt}"            
      eval ${simulation}
      # Run reconstruction 
      reconstruction="$exe_directory/ActsRecTruthTracks ${bField} --input-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_pt${pt} --output-dir reco_${detector}_e${numEvents}_t${numTracksPerEvent}_pt${pt}"
      eval ${reconstruction}
      
      echo "${jobID}, ${pt}" >> ${output_file}
      # Archive with Job ID
      mv reco_${detector}_e${numEvents}_t${numTracksPerEvent}_pt${pt}/timing.tsv timing_${jobID}.tsv
      # JobID
      let "jobID++"
   done
fi
