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
    
    if [ $suffix != truth_smeared ]; then
	    run \
  	      $outdir/performance_seeding_${suffix}.root \
    	    $refdir/performance_seeding_${suffix}.root \
      	  --title "Seeding ${suffix}" \
        	-c $config \
        	-o $outdir/seeding_${suffix}.html \
        	-p $outdir/seeding_${suffix}_plots
    fi
    
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

    Examples/Scripts/generic_plotter.py \
        $outdir/tracksummary_ckf_${suffix}.root \
        tracksummary \
        $outdir/tracksummary_ckf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/tracksummary_ckf_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/tracksummary_ckf_${suffix}.root

    run \
        $outdir/tracksummary_ckf_${suffix}_hist.root \
        $refdir/tracksummary_ckf_${suffix}_hist.root \
        --title "Track Summary CKF ${suffix}" \
        -o $outdir/tracksummary_ckf_${suffix}.html \
        -p $outdir/tracksummary_ckf_${suffix}_plots


}

full_chain truth_smeared
full_chain truth_estimated
full_chain seeded
full_chain orthogonal

run \
    $outdir/performance_gsf.root \
    $refdir/performance_gsf.root \
    --title "Truth tracking (GSF)" \
    -c CI/physmon/gsf.yml \
    -o $outdir/gsf.html \
    -p $outdir/gsf_plots

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
    $outdir/performance_ambi_orthogonal.root \
    $refdir/performance_ambi_orthogonal.root \
    --title "Ambisolver orthogonal" \
    -o $outdir/ambi_orthogonal.html \
    -p $outdir/ambi_orthogonal_plots

Examples/Scripts/vertex_mu_scan.py \
    $outdir/performance_vertexing_*mu*.root \
    $outdir/vertexing_mu_scan.pdf

rm $outdir/performance_vertexing_*mu*

exit $ec
