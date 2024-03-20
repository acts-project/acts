#!/bin/sh

bindir=$1
outdir=${2:-output}
events=${EVENTS:-32}

# ensure output directories exist
# should be handled by the tools directly
# mkdir -p ${outdir}/gen ${outdir}/sim ${outdir}/rec

# generate events with some muons
# ${bindir}/ActsGenParticleGun \
#   --events=${events} \
#   --output-dir=${outdir}/gen --output-csv=1 --output-root=1 \
#   --pg-nparticles=10 \
#   --pg-pdg=13

# simulate some events
${bindir}/ActsSimFatrasGeneric \
  --events=${events} \
  --output-dir=${outdir}/sim --output-csv=1 --output-root=1 \
  --evg-input-type=gun \
  --pg-nparticles=10 \
  --pg-pdg=13

# reconstruct simulated events
${bindir}/ActsRecTruthTracks \
  --input-dir=${outdir}/sim \
  --output-dir=${outdir}/rec --output-csv=1 --output-root=1
