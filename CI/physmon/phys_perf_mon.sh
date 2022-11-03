#!/bin/bash

set -e

outdir=$1
[ -z "$outdir" ] && outdir=physmon
mkdir -p $outdir

refdir=CI/physmon/reference
refcommit=$(cat $refdir/commit)
commit=$(git rev-parse --short HEAD)

echo "::group::Generate validation dataset"
CI/physmon/physmon.py $outdir 2>&1 > $outdir/run.log
echo "::endgroup::"

set +e

ec=0

function run() {
    a=$1
    b=$2

    echo "::group::Comparing $a vs. $b"

    histcmp \
        --label-reference=reference \
        --label-monitored=monitored \
        "$@"

    ec=$(($ec | $?))

    echo "::endgroup::"
}

function full_chain() {
    suffix=$1

    config="CI/physmon/ckf_${suffix}.yml"

    if [ ! -f "$config" ]; then
        config="CI/physmon/default.yml"
    fi
    echo $config
    
    run \
        $outdir/performance_ckf_${suffix}.root \
        $refdir/performance_ckf_${suffix}.root \
        --title "CKF ${suffix}" \
        -c $config \
        -o $outdir/ckf_${suffix}.html \
        -p $outdir/ckf_${suffix}_plots

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_vertexing_${suffix}.root \
        vertexing \
        $outdir/performance_vertexing_${suffix}_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    run \
        $outdir/performance_vertexing_${suffix}_hist.root \
        $refdir/performance_vertexing_${suffix}_hist.root \
        --title "IVF ${suffix}" \
        -o $outdir/ivf_${suffix}.html \
        -p $outdir/ivf_${suffix}_plots
}

full_chain truth_smeared
full_chain truth_estimated
full_chain seeded

run \
    $outdir/performance_truth_tracking.root \
    $refdir/performance_truth_tracking.root \
    --title "Truth tracking" \
    -c CI/physmon/truth_tracking.yml \
    -o $outdir/truth_tracking.html \
    -p $outdir/truth_tracking_plots

run \
    $outdir/performance_ambi_seeded.root \
    $refdir/performance_ambi_seeded.root \
    --title "Ambisolver seeded" \
    -o $outdir/ambi_seeded.html \
    -p $outdir/ambi_seeded_plots

run \
    $outdir/acts_analysis_residuals_and_pulls.root \
    $refdir/acts_analysis_residuals_and_pulls.root \
    --title "analysis_residuals_and_pulls" \
#    -o $outdir/analysis_residuals_and_pulls.html \
#    -p $outdir/analysis_residuals_and_pulls

Examples/Scripts/vertex_mu_scan.py \
    $outdir/performance_vertexing_*mu*.root \
    $outdir/vertexing_mu_scan.pdf

rm $outdir/performance_vertexing_*mu*

exit $ec
