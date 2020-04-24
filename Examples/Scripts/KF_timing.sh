#!/bin/bash
#
# This script runs the test of KalmanFitter timing vs. pT at different eta
# using particle gun particles
# To run the test:./KF_timing.sh -d <detector> -b <bFieldMap> -t <numTracksPerEvent> -n <numEvents>
# In default, it will run with Generic detector in constant B field using 1 event and 100 tracks/event

help()
{
   echo ""
   echo "Usage: $0 -d <detector> -b <bFieldMap> -t <numTracksPerEvent> -n <numEvents>"
   echo -e "\t-d The detector type, either 'Generic' or 'DD4hep'. Optional. In default 'Generic'" 
   echo -e "\t-b The '.txt' or '.root' file for B Field map. Optional. In default using constant BField: (0, 0, 2)"
   echo -e "\t-t The number of tracks per event. Optional. In default: 100"
   echo -e "\t-n The number of events. Optional. In default: 1"
   exit 1 # Exit script after printing help
}

#Default number of tracks/event and event
detector=Generic
numTracksPerEvent=100
numEvents=1
# Get parameters
while getopts "d:b:t:n" opt
do
   case "$opt" in
      d ) detector="$OPTARG" ;;
      b ) bFieldMap="$OPTARG" ;;
      t ) numTracksPerEvent="$OPTARG" ;;
      n ) numEvents="$OPTARG" ;;
      ? ) help ;; # Print help in case unknown parameter
   esac
done

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
echo "Test detector: ${detector}"
echo "B Field: ${bField}"
echo "Number of tracks per event: ${numTracksPerEvent}"
echo "Number of events: ${numEvents}"

time_stamp=`date +%F%T`
exe_dir=$PWD
run_dir=KalmanFitter_timing_${time_stamp}
mkdir ${run_dir}
cd ${run_dir}

output_file='output.log'
echo "*****KalmanFitter timing test vs. pT*****" > ${output_file}
echo "Test Detector: ${detector}" >> ${output_file}
echo "BField: ${bField}" >> ${output_file}
echo "Events: ${numEvents}" >> ${output_file}
echo "Tracks_per_event: ${numTracksPerEvent}" >> ${output_file} 
echo "****************************************" >> ${output_file}
echo "*"
echo "* job | eta | pt | fit_time_per_event" >> ${output_file}

jobID=0
   # Loop over the pt bins 
   for pt in 0.1 0.5 1.0 2.0 3.0 4.0 5.0 8.0 10.0 50.0 100.0 ; do
      # Loop over the eta bin number 
      for etaBin in 0 1 2 3 4; do
         etaLow=$(echo "${etaBin}*0.5"|bc) 
         etaUp=$(echo "${etaBin}*0.5 + 0.5"|bc) 
         eta=$(echo "${etaBin}*0.5 + 0.25"|bc) 

         # Run sim
         sim="${exe_dir}/ActsSimFatras${detector} ${bField} -n ${numEvents} --evg-input-type gun --pg-nparticles ${numTracksPerEvent} --pg-pt-range ${pt} ${pt} --pg-eta-range ${etaLow} ${etaUp}  --fatras-pmin-gev 0.1 --output-csv=1 --output-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}_pt${pt}"  
         echo "Run sim with '${sim}'"
         eval ${sim}
       
         # Run reco 
         reco="$exe_dir/ActsRecTruthTracks ${bField} --input-dir sim_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}_pt${pt} --output-dir reco_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}_pt${pt}"
         echo "Run reco with '${reco}'"
         eval ${reco}
         
         # Archive with Job ID
         mv reco_${detector}_e${numEvents}_t${numTracksPerEvent}_eta${eta}_pt${pt}/timing.tsv timing_${jobID}.tsv
         # Extract the fitting time
         fit_time_str=`grep "Algorithm:FittingAlgorithm" timing_${jobID}.tsv | awk '{print $3}'`
	 # Make sure the fit time is fixed-point for calculation with bc
	 fit_time_per_event=$(echo ${fit_time_str} | awk '{printf("%.10f\n", $1)}')
         fit_time_per_track=$(echo "${fit_time_per_event}/${numTracksPerEvent}"|bc -l)
         echo "${jobID}, ${etaBin}, ${pt}, ${fit_time_per_track}" >> ${output_file}
      
         # JobID
         let "jobID++"
      done
   done
