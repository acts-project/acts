#!/bin/bash

set -e


mode=${1:all}
if ! [[ $mode = @(all|kalman|gsf|fullchains|vertexing|simulation) ]]; then
    echo "Usage: $0 <all|kalman|gsf|fullchains|vertexing|simulation> (outdir)"
    exit 1
fi

outdir=${2:physmon}
[ -z "$outdir" ] && outdir=physmon
mkdir -p $outdir

refdir=CI/physmon/reference
refcommit=$(cat $refdir/commit)
commit=$(git rev-parse --short HEAD)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"  )" &> /dev/null && pwd  )

# File to accumulate the histcmp results
histcmp_results=$outdir/histcmp_results.csv

SPYRAL_BIN="spyral"
SPYRAL="${SPYRAL_BIN} run -i 0.05 --summary"

mkdir ${outdir}/memory

source $SCRIPT_DIR/setup.sh
echo "::group::Generate validation dataset"
if [[ "$mode" == "all" || "$mode" == "kalman" ]]; then
    $SPYRAL -l "Truth Tracking KF" -o "$outdir/memory/mem_truth_tracking_kalman.csv" -- CI/physmon/workflows/physmon_truth_tracking_kalman.py $outdir 2>&1 > $outdir/run_truth_tracking_kalman.log
fi
if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    $SPYRAL -l "Truth Tracking GSF" -o "$outdir/memory/mem_truth_tracking_gsf.csv" -- CI/physmon/workflows/physmon_truth_tracking_gsf.py $outdir 2>&1 > $outdir/run_truth_tracking_gsf.log
fi
if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    $SPYRAL -l "CKF Tracking" -o "$outdir/memory/mem_ckf_tracking.csv" -- CI/physmon/workflows/physmon_ckf_tracking.py $outdir 2>&1 > $outdir/run_ckf_tracking.log
fi
if [[ "$mode" == "all" || "$mode" == "vertexing" ]]; then
    $SPYRAL -l "Vertexing" -o "$outdir/memory/mem_vertexing.csv" -- CI/physmon/workflows/physmon_vertexing.py $outdir 2>&1 > $outdir/run_vertexing.log
fi
if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    $SPYRAL -l "Simulation" -o "$outdir/memory/mem_simulation.csv" -- CI/physmon/workflows/physmon_simulation.py $outdir 2>&1 > $outdir/run_simulation.log
fi
echo "::endgroup::"

$SPYRAL_BIN plot $outdir/memory/mem_truth_tracking_kalman.csv --output $outdir/memory
$SPYRAL_BIN plot $outdir/memory/mem_truth_tracking_gsf.csv --output $outdir/memory
$SPYRAL_BIN plot $outdir/memory/mem_ckf_tracking.csv --output $outdir/memory
$SPYRAL_BIN plot $outdir/memory/mem_vertexing.csv --output $outdir/memory
$SPYRAL_BIN plot $outdir/memory/mem_simulation.csv --output $outdir/memory

set +e

ec=0

function run_histcmp() {
    a=$1
    b=$2
    title=$3
    slug=$4

    echo "::group::Comparing $a vs. $b"

    if [ ! -f "$a" ]; then
      echo "::error::histcmp failed: File $a does not exist"
      ec=1
    fi

    if [ ! -f "$b" ]; then
      echo "::error::histcmp failed: File $b does not exist"
      ec=1
    fi

    histcmp \
        --label-reference=reference \
        --label-monitored=monitored \
        --title="$title" \
        -o $outdir/${slug}.html \
        -p $outdir/${slug}_plots \
        "$@"

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
      echo "::error::histcmp failed (${slug}): ec=$this_ec"
    fi

    echo "\"${title}\",${slug},${this_ec}" >> $histcmp_results

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
      run_histcmp \
        $outdir/performance_seeding_${suffix}.root \
        $refdir/performance_seeding_${suffix}.root \
        "Seeding ${suffix}" \
        -c $config \
    fi
    
    run_histcmp \
        $outdir/performance_ckf_${suffix}.root \
        $refdir/performance_ckf_${suffix}.root \
        "CKF ${suffix}" \
        ckf_${suffix} \
        -c $config

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_ivf_${suffix}.root \
        vertexing \
        $outdir/performance_ivf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_ivf_${suffix}.root

    run_histcmp \
        $outdir/performance_ivf_${suffix}_hist.root \
        $refdir/performance_ivf_${suffix}_hist.root \
        "IVF ${suffix}" \
        ivf_${suffix}

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_amvf_${suffix}.root \
        vertexing \
        $outdir/performance_amvf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_amvf_${suffix}.root

    run_histcmp \
        $outdir/performance_amvf_${suffix}_hist.root \
        $refdir/performance_amvf_${suffix}_hist.root \
        "AMVF ${suffix}" \
        amvf_${suffix}

    Examples/Scripts/generic_plotter.py \
        $outdir/tracksummary_ckf_${suffix}.root \
        tracksummary \
        $outdir/tracksummary_ckf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/tracksummary_ckf_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/tracksummary_ckf_${suffix}.root

    run_histcmp \
        $outdir/tracksummary_ckf_${suffix}_hist.root \
        $refdir/tracksummary_ckf_${suffix}_hist.root \
        "Track Summary CKF ${suffix}" \
        tracksummary_ckf_${suffix}
}

function simulation() {
    suffix=$1

    config="CI/physmon/simulation_config.yml"

    Examples/Scripts/generic_plotter.py \
        $outdir/particles_initial_${suffix}.root \
        particles \
        $outdir/particles_initial_${suffix}_hist.root \
        --silent \
        --config CI/physmon/particles_initial_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/particles_initial_${suffix}.root

    run_histcmp \
        $outdir/particles_initial_${suffix}_hist.root \
        $refdir/particles_initial_${suffix}_hist.root \
        "Particles inital ${suffix}" \
        particles_initial_${suffix}

    Examples/Scripts/generic_plotter.py \
        $outdir/particles_final_${suffix}.root \
        particles \
        $outdir/particles_final_${suffix}_hist.root \
        --silent \
        --config CI/physmon/particles_final_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/particles_final_${suffix}.root

    run_histcmp \
        $outdir/particles_final_${suffix}_hist.root \
        $refdir/particles_final_${suffix}_hist.root \
        "Particles final ${suffix}" \
        particles_final_${suffix}
}

if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    full_chain truth_smeared
    full_chain truth_estimated
    full_chain seeded
    full_chain orthogonal

    run_histcmp \
        $outdir/performance_ambi_seeded.root \
        $refdir/performance_ambi_seeded.root \
        "Ambisolver seeded"
        ambi_seeded

    run_histcmp \
        $outdir/performance_ambi_orthogonal.root \
        $refdir/performance_ambi_orthogonal.root \
        --title "Ambisolver orthogonal" \
        -o $outdir/ambi_orthogonal.html \
        -p $outdir/ambi_orthogonal_plots
fi

if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_histcmp \
        $outdir/performance_gsf.root \
        $refdir/performance_gsf.root \
        "Truth tracking (GSF)" \
        gsf \
        -c CI/physmon/gsf.yml
fi

if [[ "$mode" == "all" || "$mode" == "kalman" ]]; then
    run_histcmp \
        $outdir/performance_truth_tracking.root \
        $refdir/performance_truth_tracking.root \
        "Truth tracking" \
        truth_tracking \
        -c CI/physmon/truth_tracking.yml
fi

if [[ "$mode" == "all" || "$mode" == "vertexing" ]]; then
    Examples/Scripts/vertex_mu_scan.py \
        $outdir/performance_vertexing_*mu*.root \
        $outdir/vertexing_mu_scan.pdf

    rm $outdir/performance_vertexing_*mu*
fi

if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    simulation fatras
    simulation geant4
fi

CI/physmon/summary.py $histcmp_results \
  --base $outdir \
  --md $outdir/summary.md \
  --html $outdir/summary.html
ec=$(($ec | $?))

exit $ec
