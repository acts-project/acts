#!/bin/bash

set -e

outdir=$1
mkdir -p $outdir

refdir=CI/physmon/reference

echo "::group::Generate validation dataset"
CI/physmon/physmon.py $outdir
echo "::endgroup::"

set +e

ec=0

function run() {
    a=$1
    b=$2

    echo "::group::Comparing $a vs. $b"

    histcmp $@

    ec=$(($ec | $?))

    echo "::endgroup::"
}


run \
    physmon/performance_ckf_tracks.root \
    $refdir/performance_ckf_tracks.root \
    -o $outdir/ckf.html \


run \
    physmon/performance_truth_tracking.root \
    $refdir/performance_truth_tracking.root \
    -c CI/physmon/truth_tracking.yml \
    -o $outdir/truth_tracking.html \


exit $ec
