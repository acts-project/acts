#!/usr/bin/env python3
"""ODD KF alignment workflow (simulation / misalignment / alignment / reconstruction).

Subcommands: ``full``, ``simulation``, ``misalignment``, ``alignment``,
``reconstruction``. Shared defaults and paths are centralized in
``kfalignment_config.py`` (following the ``full_chain_odd.py`` pattern).

All sub-modules use the ``kfalignment_*.py`` naming prefix.
"""

from __future__ import annotations

import sys
from pathlib import Path

from kfalignment_paths import (
    default_digi_config,
    find_acts_source_dir,
    preload_dd4hep_for_macos,
    setup_runtime_env,
)

preload_dd4hep_for_macos()
setup_runtime_env()

import acts

from kfalignment_alignment import runAlignment
from kfalignment_config import (
    DEFAULT_NUM_EVENTS,
    build_argument_parser,
    misalignment_config_from_args,
    print_usage_tutorial,
    resolve_workflow_paths,
)
from kfalignment_misalignment import runMisalignment
from kfalignment_reconstruction import runReconstruction
from kfalignment_simulation import runSimulation

u = acts.UnitConstants
SCRIPT_DIR = Path(__file__).resolve().parent


def build_odd(misaligned: bool = False):
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector(misaligned=misaligned)
    return detector, detector.trackingGeometry()


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    if not args.command:
        print_usage_tutorial()
        sys.exit(0)

    srcdir = find_acts_source_dir()
    digi_config = default_digi_config(srcdir)
    if not digi_config.is_file():
        raise FileNotFoundError(f"Digitization config not found: {digi_config}")

    paths = resolve_workflow_paths(args, args.command, srcdir, SCRIPT_DIR, digi_config)
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    num_events = getattr(args, "num_events", DEFAULT_NUM_EVENTS)

    if args.command == "full":
        print("\n" + "=" * 80)
        print("RUNNING COMPLETE WORKFLOW")
        print("=" * 80)

        print("\n" + "=" * 80)
        print("STEP 1: SIMULATION (using NOMINAL geometry)")
        print("=" * 80)
        detector_nom, trackingGeometry_nom = build_odd(misaligned=False)
        paths.simulation_dir.mkdir(parents=True, exist_ok=True)

        runSimulation(
            trackingGeometry=trackingGeometry_nom,
            field=field,
            digiConfigFile=paths.digi_config,
            outputDir=paths.simulation_dir,
            detector=detector_nom,
            numEvents=num_events,
        ).run()

        print("\n" + "=" * 80)
        print("STEP 2a: MISALIGNMENT (write misalignment files)")
        print("=" * 80)
        paths.alignment_dir.mkdir(parents=True, exist_ok=True)
        detector_mis, trackingGeometry_mis = build_odd(misaligned=True)

        try:
            mis_cfg = misalignment_config_from_args(args)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)

        runMisalignment(
            trackingGeometry=trackingGeometry_mis,
            outputDir=paths.alignment_dir,
            config=mis_cfg,
        )

        print("\n" + "=" * 80)
        print("STEP 2b: ALIGNMENT (read misalignment; write aligned_result)")
        print("=" * 80)
        runAlignment(
            trackingGeometry=trackingGeometry_mis,
            field=field,
            inputDir=paths.simulation_dir,
            outputDir=paths.alignment_dir,
            misalignmentDir=paths.alignment_dir,
            numEvents=num_events,
            misalignmentConfig=mis_cfg,
        )

        recon_modes = ("nominal", "misaligned", "aligned")
        geo_for_mode = {
            "nominal": (detector_nom, trackingGeometry_nom),
            "misaligned": (detector_mis, trackingGeometry_mis),
            "aligned": (detector_mis, trackingGeometry_mis),
        }
        recon_dirs = {}

        for step, mode in enumerate(recon_modes, start=3):
            print("\n" + "=" * 80)
            print(f"STEP {step}: RECONSTRUCTION WITH {mode.upper()} GEOMETRY")
            print("=" * 80)
            recon_dir = paths.reconstruction_dir / mode
            recon_dir.mkdir(parents=True, exist_ok=True)
            recon_dirs[mode] = recon_dir
            _, trackingGeometry = geo_for_mode[mode]

            runReconstruction(
                trackingGeometry=trackingGeometry,
                field=field,
                inputDir=paths.simulation_dir,
                outputDir=recon_dir,
                numEvents=num_events,
                mode=mode,
                alignmentDir=paths.alignment_dir,
            ).run()

        print("\n" + "=" * 80)
        print("COMPLETE WORKFLOW FINISHED SUCCESSFULLY!")
        print("=" * 80)
        print(f"\nOutput directories:")
        print(f"  - Simulation: {paths.simulation_dir}")
        print(f"  - Alignment: {paths.alignment_dir}")
        for mode, recon_dir in recon_dirs.items():
            print(f"  - Reconstruction ({mode}): {recon_dir}")

    elif args.command == "simulation":
        paths.output_dir.mkdir(parents=True, exist_ok=True)
        detector, trackingGeometry = build_odd(misaligned=False)
        runSimulation(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=paths.digi_config,
            outputDir=paths.output_dir,
            detector=detector,
            numEvents=num_events,
        ).run()

    elif args.command == "misalignment":
        paths.output_dir.mkdir(parents=True, exist_ok=True)
        detector, trackingGeometry = build_odd(misaligned=True)
        try:
            mis_cfg = misalignment_config_from_args(args)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)

        geoIdMap, _placements = runMisalignment(
            trackingGeometry=trackingGeometry,
            outputDir=paths.output_dir,
            config=mis_cfg,
        )
        print("\n" + "=" * 60)
        print("Misalignment-only completed successfully!")
        print("=" * 60)
        print(f"\nOutput directory: {paths.output_dir}")
        print(f"Misaligned elements: {len(geoIdMap)}")
        print("Generated files:")
        print("  - misalignment_applied.txt")
        print("  - misalignment_index_map.txt")
        print(
            "\nNext: reconstruction --mode misaligned, or alignment "
            "(reads these misalignment files)"
        )

    elif args.command == "alignment":
        paths.output_dir.mkdir(parents=True, exist_ok=True)
        if not paths.input_dir.exists():
            print(f"Error: Input directory {paths.input_dir} does not exist!")
            print("Please run simulation first to generate simulation data.")
            sys.exit(1)

        detector, trackingGeometry = build_odd(misaligned=True)
        align_cfg = misalignment_config_from_args(args, target="-1")
        runAlignment(
            trackingGeometry=trackingGeometry,
            field=field,
            inputDir=paths.input_dir,
            outputDir=paths.output_dir,
            misalignmentDir=paths.misalignment_dir,
            numEvents=num_events,
            misalignmentConfig=align_cfg,
        )

    elif args.command == "reconstruction":
        paths.output_dir.mkdir(parents=True, exist_ok=True)
        if not paths.input_dir.exists():
            print(f"Error: Input directory {paths.input_dir} does not exist!")
            print("Please run simulation first to generate simulation data.")
            sys.exit(1)

        misaligned = paths.recon_mode != "nominal"
        detector, trackingGeometry = build_odd(misaligned=misaligned)

        runReconstruction(
            trackingGeometry=trackingGeometry,
            field=field,
            inputDir=paths.input_dir,
            outputDir=paths.output_dir,
            numEvents=num_events,
            mode=paths.recon_mode,
            alignmentDir=paths.alignment_dir,
        ).run()

        print("\n" + "=" * 60)
        print("Reconstruction completed successfully!")
        print("=" * 60)
        print(f"\nMode: {paths.recon_mode}")
        print(f"Output files saved to: {paths.output_dir}")
        print("\nGenerated files:")
        print("  - trackstates_kf.root: Track states")
        print("  - tracksummary_kf.root: Track summary")
        print("  - performance_kf.root: Reconstruction performance")


if __name__ == "__main__":
    main()
