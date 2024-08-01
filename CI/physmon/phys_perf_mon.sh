#!/bin/bash

set -e

# helper function to selectively print and run commands without a subshell
function run() {
    set -x
    "$@"
    # save exit code
    { rec=$?; } 2> /dev/null
    { set +x;   } 2> /dev/null
    # restore exit code
    (exit $rec)
}

export run


shopt -s extglob


mode=${1:-all}
if ! [[ $mode = @(all|kf|gsf|gx2f|fullchains|simulation) ]]; then
    echo "Usage: $0 <all|kf|gsf|gx2f|fullchains|simulation> (outdir)"
    exit 1
fi

outdir=${2:-physmon}
mkdir -p $outdir

refdir=CI/physmon/reference
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"  )" &> /dev/null && pwd  )

# File to accumulate the histcmp results
histcmp_results=$outdir/histcmp_results.csv
echo -n "" > $histcmp_results

memory_dir=${outdir}/memory
mkdir -p "$memory_dir"

if [ "$(uname)" == "Darwin" ]; then
    function measure {
        label=$1
        slug=$2
        shift
        shift
        echo "Measure Darwin $label ($slug)" >&2
        tmp=$(mktemp)
        echo "+ $@" >&2
        /usr/bin/time -l -o "$tmp" "$@"
        # save exit code
        rec=$?

        of="${memory_dir}/mem_${slug}.csv"
        {
            echo "# spyral-label: $label"
            echo "# spyral-cmd: $*"
            echo "time,rss,vms"
            # in bytes
            grep -E "real" "$tmp" | awk '{printf $1}'
            printf ","
            grep -E "maximum resident set size" "$tmp" | awk '{printf $1}'
            printf ",0\n"
        } > "$of"
        # restore exit code
        (exit $rec)
    }
    export measure
elif [ "$(uname)" == "Linux" ]; then
    function measure {
        label=$1
        slug=$2
        shift
        shift
        echo "Measure Linux $label ($slug)" >&2
        tmp=$(mktemp)
        echo "+ $@" >&2
        /usr/bin/time -v -o  "$tmp" "$@"
        # save exit code
        rec=$?
        # in kbytes
        max_rss=$(grep "Maximum resident set size (kbytes):" "$tmp" | awk '{printf $(NF)}')
        max_rss=$(( 1000*max_rss ))
        wall_time=$(grep "Elapsed (wall clock)" "$tmp" | awk '{printf $(NF)}')
        echo $max_rss
        wall_time=$(python3 -c "i='${wall_time}';p=i.split(':');p = p if len(p) == 3 else ['0', *p];t=float(p[0])*60*60 + float(p[1])*60 + float(p[2]);print(t)")
        echo $wall_time

        of="${memory_dir}/mem_${slug}.csv"
        {
            echo "# spyral-label: $label"
            echo "# spyral-cmd: $*"
            echo "time,rss,vms"
            echo "${wall_time},${max_rss},0"
        } > "$of"
        # restore exit code
        (exit $rec)
    }
    export measure
else
    function measure {
        echo "Not measuring because unknown environment"
        shift
        shift
        "$@"
    }
    export measure
fi

if [ -n "$CI" ]; then
    echo "CI mode, do not abort immediately on failure"
    set +e
fi
ec=0

source $SCRIPT_DIR/setup.sh

function run_physmon_gen() {
    title=$1
    slug=$2

    script=CI/physmon/workflows/physmon_${slug}.py

    mkdir -p $outdir/data/$slug
    measure "$title" "$slug" ${script} $outdir/data/$slug 2>&1 > $outdir/data/$slug/run_${slug}.log

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
        echo "::error::🟥 Dataset generation failed: ${script} -> ec=$this_ec"
    else
        echo "::notice::✅ Dataset generation succeeded: ${script}"
    fi
}

echo "::group::Generate validation dataset"
if [[ "$mode" == "all" || "$mode" == "kf" ]]; then
    run_physmon_gen "Truth Tracking KF" "trackfitting_kf"
fi
if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_physmon_gen "Truth Tracking GSF" "trackfitting_gsf"
fi
if [[ "$mode" == "all" || "$mode" == "gx2f" ]]; then
    run_physmon_gen "Truth Tracking GX2F" "trackfitting_gx2f"
fi
if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    run_physmon_gen "CKF single muon" "trackfinding_1muon"
    run_physmon_gen "CKF muon 50" "trackfinding_4muon_50vertices"
    run_physmon_gen "CKF ttbar 200" "trackfinding_ttbar_pu200"
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

    run histcmp $a $b \
        --label-reference=reference \
        --label-monitored=monitored \
        --title="$title" \
        -o $outdir/html/${slug}.html \
        -p $outdir/html/${slug}_plots \
        "$@"

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
        echo "::error::histcmp failed (${slug}): ec=$this_ec"
    fi

    echo "\"${title}\",${slug},${this_ec}" >> $histcmp_results

    echo "::endgroup::"
}

function trackfinding() {
    slug=$1

    config="CI/physmon/config/default.yml"

    if [ -f $refdir/$slug/performance_seeding.root ]; then
        run_histcmp \
            $outdir/data/$slug/performance_seeding.root \
            $refdir/$slug/performance_seeding.root \
            "Seeding ${slug}" \
            seeding_${slug} \
            -c $config
    fi

    run_histcmp \
        $outdir/data/$slug/performance_ckf.root \
        $refdir/$slug/performance_ckf.root \
        "CKF ${slug}" \
        ckf_${slug} \
        -c $config

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$slug/tracksummary_ckf.root \
        tracksummary \
        $outdir/data/$slug/tracksummary_ckf_hist.root \
        --silent \
        --config CI/physmon/config/tracksummary_ckf.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/$slug/tracksummary_ckf.root

    run_histcmp \
        $outdir/data/$slug/tracksummary_ckf_hist.root \
        $refdir/$slug/tracksummary_ckf_hist.root \
        "Track Summary CKF ${slug}" \
        tracksummary_ckf_${slug}

    if [ -f $refdir/$slug/performance_ambi.root ]; then
        run_histcmp \
            $outdir/data/$slug/performance_ambi.root \
            $refdir/$slug/performance_ambi.root \
            "Ambisolver ${slug}" \
            ambi_${slug}
    fi
}

function vertexing() {
    slug=$1
    config=$2

    if [ -f $refdir/$slug/performance_ivf_notime_hist.root ]; then
        run Examples/Scripts/generic_plotter.py \
            $outdir/data/$slug/performance_ivf_notime.root \
            vertexing \
            $outdir/data/$slug/performance_ivf_notime_hist.root \
            --silent \
            --config $config
        ec=$(($ec | $?))

        run_histcmp \
            $outdir/data/$slug/performance_ivf_notime_hist.root \
            $refdir/$slug/performance_ivf_notime_hist.root \
            "IVF notime ${slug}" \
            ivf_notime_${slug}
    fi

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$slug/performance_amvf_gauss_notime.root \
        vertexing \
        $outdir/data/$slug/performance_amvf_gauss_notime_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    run_histcmp \
        $outdir/data/$slug/performance_amvf_gauss_notime_hist.root \
        $refdir/$slug/performance_amvf_gauss_notime_hist.root \
        "AMVF gauss notime ${slug}" \
        amvf_gauss_notime_${slug}

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$slug/performance_amvf_grid_time.root \
        vertexing \
        $outdir/data/$slug/performance_amvf_grid_time_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    run_histcmp \
        $outdir/data/$slug/performance_amvf_grid_time_hist.root \
        $refdir/$slug/performance_amvf_grid_time_hist.root \
        "AMVF grid time ${slug}" \
        amvf_grid_time_${slug}
}

function simulation() {
    suffix=$1

    config="CI/physmon/config/simulation.yml"

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/particles_${suffix}.root \
        particles \
        $outdir/data/simulation/particles_${suffix}_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/simulation/particles_${suffix}.root

    run_histcmp \
        $outdir/data/simulation/particles_${suffix}_hist.root \
        $refdir/simulation/particles_${suffix}_hist.root \
        "Particles ${suffix}" \
        particles_${suffix}
}

function generation() {
    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/pythia8_particles_ttbar.root \
        particles \
        $outdir/data/simulation/particles_ttbar_hist.root \
        --silent \
        --config CI/physmon/config/pythia8_ttbar.yml

    run_histcmp \
        $outdir/data/simulation/particles_ttbar_hist.root \
        $refdir/simulation/particles_ttbar_hist.root \
        "Particles ttbar" \
        particles_ttbar

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/pythia8_vertices_ttbar.root \
        vertices \
        $outdir/data/simulation/vertices_ttbar_hist.root \
        --silent \
        --config CI/physmon/config/pythia8_ttbar.yml

    run_histcmp \
        $outdir/data/simulation/vertices_ttbar_hist.root \
        $refdir/simulation/vertices_ttbar_hist.root \
        "Vertices ttbar" \
        vertices_ttbar
}

if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    simulation fatras
    simulation geant4

    generation
fi

if [[ "$mode" == "all" || "$mode" == "kf" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_kf/performance_kf.root \
        $refdir/trackfitting_kf/performance_kf.root \
        "Truth tracking (KF)" \
        kf \
        -c CI/physmon/config/trackfitting_kf.yml
fi

if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_gsf/performance_gsf.root \
        $refdir/trackfitting_gsf/performance_gsf.root \
        "Truth tracking (GSF)" \
        gsf \
        -c CI/physmon/config/trackfitting_gsf.yml
fi

if [[ "$mode" == "all" || "$mode" == "gx2f" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_gx2f/performance_gx2f.root \
        $refdir/trackfitting_gx2f/performance_gx2f.root \
        "Truth tracking (GX2F)" \
        gx2f \
        -c CI/physmon/config/trackfitting_gx2f.yml
fi

if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    trackfinding trackfinding_1muon/truth_smeared
    trackfinding trackfinding_1muon/truth_estimated
    trackfinding trackfinding_1muon/seeded
    trackfinding trackfinding_1muon/orthogonal

    trackfinding trackfinding_4muon_50vertices
    vertexing trackfinding_4muon_50vertices CI/physmon/config/vertexing_4muon_50vertices.yml

    trackfinding trackfinding_ttbar_pu200
    vertexing trackfinding_ttbar_pu200 CI/physmon/config/vertexing_ttbar_pu200.yml
fi

run CI/physmon/summary.py $histcmp_results \
  --md $outdir/summary.md \
  --html $outdir/summary.html
ec=$(($ec | $?))

exit $ec
