#!/bin/bash

set -e

outdir=$1
mkdir -p $outdir

# GENERATE THE OUTPUTS!

histcmp \
    histcmp_demo/performance_track_fitter_a.root \
    histcmp_demo/performance_track_fitter_b.root \
    -o $outdir/track_fitter.html \

histcmp \
    histcmp_demo/performance_ckf.root \
    histcmp_demo/performance_ckf_main.root \
    -o $outdir/ckf.html \

