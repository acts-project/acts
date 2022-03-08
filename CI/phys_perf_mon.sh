#!/bin/bash

set -e

outdir=$1
mkdir -p $outdir

# GENERATE THE OUTPUTS!

set +e

ec=0

histcmp \
    histcmp_demo/performance_ckf.root \
    histcmp_demo/performance_ckf_main.root \
    -o $outdir/ckf.html \

ec=$(($ec | $?))

histcmp \
    histcmp_demo/performance_track_fitter_a.root \
    histcmp_demo/performance_track_fitter_b.root \
    -c track_fitter.yml \
    -o $outdir/track_fitter.html \

ec=$(($ec | $?))


exit $ec
