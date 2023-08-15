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

SPYRAL_BIN="spyral"
SPYRAL="${SPYRAL_BIN} run -i 0.1 --summary"

mkdir -p "${outdir}/memory"

source $SCRIPT_DIR/setup.sh
echo "::group::Generate validation dataset"
if [[ "$mode" == "all" || "$mode" == "kalman" ]]; then
    $SPYRAL -l "Truth Tracking KF" -o "$outdir/memory/mem_truth_tracking_kalman.csv" -- CI/physmon/workflows/physmon_truth_tracking_kalman.py $outdir 2>&1 > $outdir/run_truth_tracking_kalman.log
    $SPYRAL_BIN plot $outdir/memory/mem_truth_tracking_kalman.csv --output $outdir/memory
fi
if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    $SPYRAL -l "Truth Tracking GSF" -o "$outdir/memory/mem_truth_tracking_gsf.csv" -- CI/physmon/workflows/physmon_truth_tracking_gsf.py $outdir 2>&1 > $outdir/run_truth_tracking_gsf.log
    $SPYRAL_BIN plot $outdir/memory/mem_truth_tracking_gsf.csv --output $outdir/memory
fi
if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    $SPYRAL -l "CKF Tracking" -o "$outdir/memory/mem_ckf_tracking.csv" -- CI/physmon/workflows/physmon_ckf_tracking.py $outdir 2>&1 > $outdir/run_ckf_tracking.log
    $SPYRAL -l "Track finding ttbar" -o "$outdir/memory/mem_ttbar.csv" -- CI/physmon/workflows/physmon_track_finding_ttbar.py $outdir 2>&1 > $outdir/run_track_finding_ttbar.log
    $SPYRAL_BIN plot $outdir/memory/mem_ckf_tracking.csv --output $outdir/memory
    $SPYRAL_BIN plot $outdir/memory/mem_ttbar.csv --output $outdir/memory
fi
if [[ "$mode" == "all" || "$mode" == "vertexing" ]]; then
    $SPYRAL -l "Vertexing" -o "$outdir/memory/mem_vertexing.csv" -- CI/physmon/workflows/physmon_vertexing.py $outdir 2>&1 > $outdir/run_vertexing.log
    $SPYRAL_BIN plot $outdir/memory/mem_vertexing.csv --output $outdir/memory
fi
if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    $SPYRAL -l "Simulation" -o "$outdir/memory/mem_simulation.csv" -- CI/physmon/workflows/physmon_simulation.py $outdir 2>&1 > $outdir/run_simulation.log
    $SPYRAL_BIN plot $outdir/memory/mem_simulation.csv --output $outdir/memory
fi
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
        $outdir/performance_ivf_${suffix}.root \
        vertexing \
        $outdir/performance_ivf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_ivf_${suffix}.root

    run \
        $outdir/performance_ivf_${suffix}_hist.root \
        $refdir/performance_ivf_${suffix}_hist.root \
        --title "IVF ${suffix}" \
        -o $outdir/ivf_${suffix}.html \
        -p $outdir/ivf_${suffix}_plots

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_amvf_${suffix}.root \
        vertexing \
        $outdir/performance_amvf_${suffix}_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_amvf_${suffix}.root

    run \
        $outdir/performance_amvf_${suffix}_hist.root \
        $refdir/performance_amvf_${suffix}_hist.root \
        --title "AMVF ${suffix}" \
        -o $outdir/amvf_${suffix}.html \
        -p $outdir/amvf_${suffix}_plots

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

    run \
        $outdir/particles_initial_${suffix}_hist.root \
        $refdir/particles_initial_${suffix}_hist.root \
        --title "Particles inital ${suffix}" \
        -o $outdir/particles_initial_${suffix}.html \
        -p $outdir/particles_initial_${suffix}_plots

    Examples/Scripts/generic_plotter.py \
        $outdir/particles_final_${suffix}.root \
        particles \
        $outdir/particles_final_${suffix}_hist.root \
        --silent \
        --config CI/physmon/particles_final_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/particles_final_${suffix}.root

    run \
        $outdir/particles_final_${suffix}_hist.root \
        $refdir/particles_final_${suffix}_hist.root \
        --title "Particles final ${suffix}" \
        -o $outdir/particles_final_${suffix}.html \
        -p $outdir/particles_final_${suffix}_plots
}

if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    full_chain truth_smeared
    full_chain truth_estimated
    full_chain seeded
    full_chain orthogonal

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

    run \
        $outdir/performance_seeding_ttbar.root \
        $refdir/performance_seeding_ttbar.root \
        --title "Seeding ttbar" \
        -c $config \
        -o $outdir/seeding_ttbar.html \
        -p $outdir/seeding_ttbar_plots

    run \
        $outdir/performance_ckf_ttbar.root \
        $refdir/performance_ckf_ttbar.root \
        --title "CKF ttbar" \
        -c $config \
        -o $outdir/ckf_ttbar.html \
        -p $outdir/ckf_ttbar_plots

    run \
        $outdir/performance_ambi_ttbar.root \
        $refdir/performance_ambi_ttbar.root \
        --title "Ambisolver " \
        -o $outdir/ambi_ttbar.html \
        -p $outdir/ambi_ttbar_plots

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_amvf_ttbar.root \
        vertexing \
        $outdir/performance_amvf_ttbar_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_amvf_ttbar.root

    run \
        $outdir/performance_amvf_ttbar_hist.root \
        $refdir/performance_amvf_ttbar_hist.root \
        --title "AMVF ttbar" \
        -o $outdir/amvf_ttbar.html \
        -p $outdir/amvf_ttbar_plots
fi

if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run \
        $outdir/performance_gsf.root \
        $refdir/performance_gsf.root \
        --title "Truth tracking (GSF)" \
        -c CI/physmon/gsf.yml \
        -o $outdir/gsf.html \
        -p $outdir/gsf_plots
fi

if [[ "$mode" == "all" || "$mode" == "kalman" ]]; then
    run \
        $outdir/performance_truth_tracking.root \
        $refdir/performance_truth_tracking.root \
        --title "Truth tracking" \
        -c CI/physmon/truth_tracking.yml \
        -o $outdir/truth_tracking.html \
        -p $outdir/truth_tracking_plots
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

CI/physmon/summary.py $outdir/*.html \
  --base $outdir \
  --md $outdir/summary.md \
  --html $outdir/summary.html
ec=$(($ec | $?))

exit $ec
