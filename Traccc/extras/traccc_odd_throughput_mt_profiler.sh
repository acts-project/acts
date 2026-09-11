#!/bin/bash
#
# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Simple script running the selected instance of the multi-threaded throughput
# executable on a whole set of ODD ttbar simulations, with different pileup
# values.
#

# Stop on errors.
set -e

# Function printing the usage information for the script.
usage() {
   echo "Script running a suite of multi-threaded throughput tests"
   echo ""
   echo "Usage: traccc_throughput_mt_profiling.sh [options]"
   echo ""
   echo "Basic options:"
   echo "  -x <executable>      Selects the executable to use"
   echo "  -i <inputDir>        Selects the input directory with the ttbar"
   echo "                       simulations"
   echo "  -m <minThreads>      Minimum number of threads to test"
   echo "  -t <maxThreads>      Maximum number of threads to test"
   echo "  -s <threadStep>      Steps to increase the thread count by"
   echo "  -r <repetitions>     The number of repetitions in the test"
   echo "  -e <eventMultiplier> Multiplier for the number of events per thread"
   echo "  -c <csvFile>         Name of the output CSV file"
   echo "  -h                   Print this help"
   echo ""
}

# Parse the command line arguments.
TRACCC_EXECUTABLE=${TRACCC_EXECUTABLE:-"traccc_throughput_mt"}
TRACCC_INPUT_DIR=${TRACCC_INPUT_DIR:-"odd/"}
TRACCC_MIN_THREADS=${TRACCC_MIN_THREADS:-1}
TRACCC_MAX_THREADS=${TRACCC_MAX_THREADS:-$(nproc)}
TRACCC_THREAD_STEP=${TRACCC_THREAD_STEP:-1}
TRACCC_REPETITIONS=${TRACCC_REPETITIONS:-5}
TRACCC_EVT_MULTI=${TRACCC_EVT_MULTI:-"1"}
TRACCC_CSV_FILE=${TRACCC_CSV_FILE:-"output.csv"}
while getopts ":x:i:m:t:r:e:c:h" opt; do
   case $opt in
      x)
         TRACCC_EXECUTABLE=$OPTARG
         ;;
      i)
         TRACCC_INPUT_DIR=$OPTARG
         ;;
      m)
         TRACCC_MIN_THREADS=$OPTARG
         ;;
      t)
         TRACCC_MAX_THREADS=$OPTARG
         ;;
      s)
         TRACCC_THREAD_STEP=$OPTARG
         ;;
      r)
         TRACCC_REPETITIONS=$OPTARG
         ;;
      e)
         TRACCC_EVT_MULTI=$OPTARG
         ;;
      c)
         TRACCC_CSV_FILE=$OPTARG
         ;;
      h)
         usage
         exit 0
         ;;
      :)
         echo "Argument -$OPTARG requires a parameter!"
         usage
         exit 1
         ;;
      ?)
         echo "Unknown argument: -$OPTARG"
         usage
         exit 1
         ;;
   esac
done

# Print the configuration received.
echo "Using configuration:"
echo "   EXECUTABLE  : ${TRACCC_EXECUTABLE}"
echo "   INPUT_DIR   : ${TRACCC_INPUT_DIR}"
echo "   MIN_THREADS : ${TRACCC_MIN_THREADS}"
echo "   MAX_THREADS : ${TRACCC_MAX_THREADS}"
echo "   THREAD_STEP : ${TRACCC_THREAD_STEP}"
echo "   REPETITIONS : ${TRACCC_REPETITIONS}"
echo "   EVT_MULTI   : ${TRACCC_EVT_MULTI}"
echo "   CSV_FILE    : ${TRACCC_CSV_FILE}"

# Check whether the output file already exists. Refuse to overwrite existing
# files.
if [[ -f "${TRACCC_CSV_FILE}" ]]; then
   echo "***"
   echo "*** Will not overwrite ${TRACCC_CSV_FILE}!"
   echo "***"
   exit 1
fi

# The input directories to use.
TRACCC_INPUT_DIRS=("ttbar_mu20"  "ttbar_mu40"  "ttbar_mu60"  "ttbar_mu80"
                   "ttbar_mu100" "ttbar_mu140" "ttbar_mu200" "ttbar_mu300")

# The number of events to process for the different mu values. Chosen to take
# roughly the same amount of time to process on a CPU.
declare -A TRACCC_EVT_COUNT
TRACCC_EVT_COUNT["ttbar_mu20"]=$((50*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu40"]=$((50*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu60"]=$((25*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu80"]=$((25*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu100"]=$((20*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu140"]=$((15*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu200"]=$((8*${TRACCC_EVT_MULTI}))
TRACCC_EVT_COUNT["ttbar_mu300"]=$((5*${TRACCC_EVT_MULTI}))

# Put a header on the CSV file.
echo "directory,threads,loaded_events,cold_run_events,processed_events,warm_up_time,processing_time" \
   > "${TRACCC_CSV_FILE}"

# Counter for a nice printout.
COUNTER=1
COUNT=$((${#TRACCC_INPUT_DIRS[@]}*${TRACCC_MAX_THREADS}*${TRACCC_REPETITIONS}))

# Iterate over the number of threads.
for NTHREAD in $(seq ${TRACCC_MIN_THREADS} ${TRACCC_THREAD_STEP} ${TRACCC_MAX_THREADS}); do
   # Iterate over the input datasets.
   for EVTDIR in ${TRACCC_INPUT_DIRS[@]}; do
      # Perform the requested number of repetitions.
      for REPEAT in $(seq ${TRACCC_REPETITIONS}); do

         # Tell the user what's happening.
         echo ""
         echo "Running test ${COUNTER} / ${COUNT}"
         ((COUNTER++))

         # Run the throughput test.
         ${TRACCC_EXECUTABLE}                                                   \
            --detector-file="geometries/odd/odd-detray_geometry_detray.json"    \
            --material-file="geometries/odd/odd-detray_material_detray.json"    \
            --grid-file="geometries/odd/odd-detray_surface_grids_detray.json"   \
            --digitization-file="geometries/odd/odd-digi-geometric-config.json" \
            --conditions-file="geometries/odd/odd-conditions.json"               \
            --input-directory="${TRACCC_INPUT_DIR}/geant4_${EVTDIR}/"           \
            --input-events=500                                                  \
            --cpu-threads=${NTHREAD}                                            \
            --cold-run-events=$((5*${NTHREAD}))                                 \
            --processed-events=$((${TRACCC_EVT_COUNT[${EVTDIR}]}*${NTHREAD}))   \
            --log-file="${TRACCC_CSV_FILE}"
      done
   done
done
