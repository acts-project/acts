#!/bin/bash

set -e
set -x


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
echo -n "" > $histcmp_results

SPYRAL_BIN=${SPYRAL_BIN:=spyral}
SPYRAL="${SPYRAL_BIN} run -i 0.05 --summary"

mkdir -p ${outdir}/memory

set +e
ec=0

source $SCRIPT_DIR/setup.sh

function run_physmon_gen() {
    title=$1
    slug=$2

    script=CI/physmon/workflows/physmon_${slug}.py

    $SPYRAL -l "$title" -o "$outdir/memory/mem_${slug}.csv" -- ${script} $outdir 2>&1 > $outdir/run_${slug}.log

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
      echo "::error::🟥 Dataset generation failed: ${script} -> ec=$this_ec"
    else
      echo "::notice::✅ Dataset generation succeeded: ${script}"   
    fi

    $SPYRAL_BIN plot $outdir/memory/mem_${slug}.csv --output $outdir/memory
}

echo "::group::Generate validation dataset"
if [[ "$mode" == "all" || "$mode" == "kalman" ]]; then
    run_physmon_gen "Truth Tracking KF" "truth_tracking_kalman"
fi
if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_physmon_gen "Truth Tracking GSF" "truth_tracking_gsf"
fi
if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    run_physmon_gen "CKF Tracking" "ckf_tracking"
    run_physmon_gen "Track finding ttbar" "track_finding_ttbar"

fi
if [[ "$mode" == "all" || "$mode" == "vertexing" ]]; then
    run_physmon_gen "Vertexing" "vertexing"
fi
if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    run_physmon_gen "Simulation" "simulation"
fi
echo "::endgroup::"


function run_histcmp() {
    a=$1
    b=$2
    title=$3
    slug=$4
    shift 4

    echo "::group::Comparing $a vs. $b"

    if [ ! -f "$a" ]; then
      echo "::error::histcmp failed: File $a does not exist"
      ec=1
    fi

    if [ ! -f "$b" ]; then
      echo "::error::histcmp failed: File $b does not exist"
      ec=1
    fi

    histcmp $a $b \
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
        seeding_${suffix} \
        -c $config
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

    if [ $suffix == seeded ]; then
	    Examples/Scripts/generic_plotter.py \
            $outdir/performance_amvf_gridseeder_${suffix}.root \
            vertexing \
            $outdir/performance_amvf_gridseeder_${suffix}_hist.root \
            --silent \
            --config CI/physmon/vertexing_config.yml
        ec=$(($ec | $?))

        # remove ntuple file because it's large
        rm $outdir/performance_amvf_gridseeder_${suffix}.root

        run_histcmp \
            $outdir/performance_amvf_gridseeder_${suffix}_hist.root \
            $refdir/performance_amvf_gridseeder_${suffix}_hist.root \
            "AMVF (+grid seeder) ${suffix}" \
            amvf_gridseeder_${suffix}
    fi

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
        "Ambisolver seeded" \
        ambi_seeded

    run_histcmp \
        $outdir/performance_ambi_orthogonal.root \
        $refdir/performance_ambi_orthogonal.root \
        "Ambisolver orthogonal" \
        ambi_orthogonal

    run_histcmp \
        $outdir/performance_seeding_ttbar.root \
        $refdir/performance_seeding_ttbar.root \
        "Seeding ttbar" \
        seeding_ttbar \
        -c $config

    run_histcmp \
        $outdir/performance_ckf_ttbar.root \
        $refdir/performance_ckf_ttbar.root \
        "CKF ttbar" \
        ckf_ttbar \
        -c $config

    run_histcmp \
        $outdir/performance_ambi_ttbar.root \
        $refdir/performance_ambi_ttbar.root \
        "Ambisolver " \
        ambi_ttbar

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_amvf_ttbar.root \
        vertexing \
        $outdir/performance_amvf_ttbar_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    Examples/Scripts/generic_plotter.py \
        $outdir/tracksummary_ckf_ttbar.root \
        tracksummary \
        $outdir/tracksummary_ckf_ttbar_hist.root \
        --config CI/physmon/tracksummary_ckf_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/tracksummary_ckf_ttbar.root

    run_histcmp \
        $outdir/tracksummary_ckf_ttbar_hist.root \
        $refdir/tracksummary_ckf_ttbar_hist.root \
        "Track Summary CKF ttbar" \
        tracksummary_ckf_ttbar

    # remove ntuple file because it's large
    rm $outdir/performance_amvf_ttbar.root

    run_histcmp \
        $outdir/performance_amvf_ttbar_hist.root \
        $refdir/performance_amvf_ttbar_hist.root \
        "AMVF ttbar" \
        amvf_ttbar

    Examples/Scripts/generic_plotter.py \
        $outdir/performance_amvf_gridseeder_ttbar.root \
        vertexing \
        $outdir/performance_amvf_gridseeder_ttbar_hist.root \
        --silent \
        --config CI/physmon/vertexing_config.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/performance_amvf_gridseeder_ttbar.root

    run_histcmp \
        $outdir/performance_amvf_gridseeder_ttbar_hist.root \
        $refdir/performance_amvf_gridseeder_ttbar_hist.root \
        "AMVF (+grid seeder) ttbar" \
        amvf_gridseeder_ttbar
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
  --md $outdir/summary.md \
  --html $outdir/summary.html
ec=$(($ec | $?))

exit $ec
