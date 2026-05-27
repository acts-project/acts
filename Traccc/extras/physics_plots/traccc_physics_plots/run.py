import subprocess
import logging
import re
import shutil
import os
import json
import pathlib

log = logging.getLogger("traccc_physics_plots.run")

SEEDING_EXAMPLE_ARGS = [
    "--input-directory=/data/Acts/odd-simulations-20240506/geant4_ttbar_mu200",
    "--digitization-file=geometries/odd/odd-digi-geometric-config.json",
    "--conditions-file=geometries/odd/odd-conditions.json",
    "--detector-file=geometries/odd/odd-detray_geometry_detray.json",
    "--grid-file=geometries/odd/odd-detray_surface_grids_detray.json",
    "--material-file=geometries/odd/odd-detray_material_detray.json",
    "--input-events=10",
    "--use-acts-geom-source=on",
    "--check-performance",
    "--truth-finding-min-track-candidates=5",
    "--truth-finding-min-pt=1.0",
    "--truth-finding-min-z=-150",
    "--truth-finding-max-z=150",
    "--truth-finding-max-r=10",
    "--seed-matching-ratio=0.99",
    "--track-matching-ratio=0.5",
    "--track-candidates-range=5:100",
    "--seedfinder-vertex-range=-150:150",
]


def run_convert_seeding_example(
    seeding_example_executable,
    seeding_example_args,
    tmpdirname,
    output_dir,
    root_to_csv_exec=None,
):
    log.info("Running seeding example...")

    if root_to_csv_exec is None:
        root_to_csv_exec = "./root_to_csv/root_to_csv"

    if seeding_example_executable.name != "traccc_seeding_example_cuda":
        log.warning(
            'Executable name is not "traccc_seeding_example_cuda" but "%s"; this will likely fail',
            seeding_example_executable.name,
        )

    result = subprocess.run(
        [seeding_example_executable] + seeding_example_args,
        check=True,
        cwd=tmpdirname,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    stdout = result.stdout.decode("utf-8")

    if match := re.search(r"- created \(cuda\)\s+(\d+) seeds", stdout):
        num_seeds = int(match.group(1))
    else:
        raise ValueError("Seed count could not be parsed from stdout!")

    if match := re.search(r"- created \(cuda\)\s+(\d+) found tracks", stdout):
        num_found_tracks = int(match.group(1))
    else:
        raise ValueError("Found track count could not be parsed from stdout!")

    if match := re.search(r"- created \(cuda\)\s+(\d+) fitted tracks", stdout):
        num_fitted_tracks = int(match.group(1))
    else:
        raise ValueError("Fitted track count could not be parsed from stdout!")

    with open(output_dir / "counts.json", "w") as f:
        json.dump(
            {
                "seeds": num_seeds,
                "found": num_found_tracks,
                "fitted": num_fitted_tracks,
            },
            f,
        )

    for f in ["seeding", "finding", "fitting"]:
        fn = "performance_track_%s.root" % f
        shutil.copy(pathlib.Path(tmpdirname) / fn, output_dir)

        try:
            os.mkdir(output_dir / f)
        except FileExistsError:
            pass

        subprocess.run(
            [
                root_to_csv_exec,
                str(output_dir / fn),
                str(output_dir / f),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
