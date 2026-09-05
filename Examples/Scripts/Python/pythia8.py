#!/usr/bin/env python3
import argparse
from pathlib import Path

import acts
import acts.examples
import acts.examples.hepmc3
from acts.examples.simulation import addPythia8

u = acts.UnitConstants

DEFAULT_HARD_PROCESS = ["Top:qqbar2ttbar=on"]
DEFAULT_HEPMC3_FILENAME = "events.hepmc3"
HEPMC3_COMPRESSION_CHOICES = tuple(
    c.name for c in acts.examples.hepmc3.availableCompressionModes()
)


def runPythia8(
    outputDir,
    outputRoot: bool = True,
    outputCsv: bool = False,
    s: acts.examples.Sequencer = None,
    vtxGen=None,
    hardProcess=None,
    outputHepMC3: bool = False,
    hepmc3Compression=acts.examples.hepmc3.Compression.none,
):
    # Preliminaries
    rnd = acts.examples.RandomNumbers()
    outputDir = Path(outputDir)

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    addPythia8(
        s,
        rnd=rnd,
        hardProcess=hardProcess,
        outputDirCsv=outputDir / "csv" if outputCsv else None,
        outputDirRoot=outputDir if outputRoot else None,
        writeHepMC3=None,
        vtxGen=vtxGen,
    )

    if outputHepMC3:
        outputHepMC3Path = outputDir / DEFAULT_HEPMC3_FILENAME
        s.addWriter(
            acts.examples.hepmc3.HepMC3Writer(
                acts.logging.INFO,
                inputEvent="pythia8-event",
                outputPath=outputHepMC3Path,
                compression=hepmc3Compression,
            )
        )

    return s


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Run Pythia8 event generation with optional HepMC3 output."
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="Output directory for ROOT outputs.",
    )
    parser.add_argument(
        "-n",
        "--events",
        type=int,
        default=10,
        help="Number of events to run.",
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=-1,
        help="Number of sequencer threads (-1 for auto).",
    )
    parser.add_argument(
        "--no-root",
        action="store_true",
        help="Disable ROOT output.",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Write CSV output to output-dir/csv.",
    )
    parser.add_argument(
        "--hepmc3",
        action="store_true",
        help=f"Write HepMC3 using HepMC3Writer to output-dir/{DEFAULT_HEPMC3_FILENAME}.",
    )
    parser.add_argument(
        "--hepmc3-compression",
        choices=HEPMC3_COMPRESSION_CHOICES,
        default=("zstd" if "zstd" in HEPMC3_COMPRESSION_CHOICES else "none"),
        help="Compression mode for HepMC3 output.",
    )
    parser.add_argument(
        "-p",
        "--process",
        action="append",
        default=None,
        help=(
            "Pythia process setting; can be provided multiple times. "
            'Defaults to ttbar: "Top:qqbar2ttbar=on".'
        ),
    )
    parser.add_argument(
        "--beamspot-mean-mm",
        nargs=3,
        type=float,
        default=(0.0, 0.0, 0.0),
        metavar=("X", "Y", "Z"),
        help="Beamspot mean position in mm.",
    )
    parser.add_argument(
        "--beamspot-sigma-mm",
        nargs=3,
        type=float,
        default=(0.0, 0.0, 0.0),
        metavar=("SX", "SY", "SZ"),
        help="Beamspot Gaussian sigma in mm.",
    )
    parser.add_argument(
        "--beamspot-mean-time-ns",
        type=float,
        default=0.0,
        help="Beamspot mean time in ns.",
    )
    parser.add_argument(
        "--beamspot-sigma-time-ns",
        type=float,
        default=0.0,
        help="Beamspot Gaussian sigma in ns.",
    )
    return parser.parse_args()


if "__main__" == __name__:
    args = _parse_args()
    outputDir = args.output_dir
    outputDir.mkdir(parents=True, exist_ok=True)
    hardProcess = args.process if args.process is not None else DEFAULT_HARD_PROCESS
    hepmc3Compression = getattr(
        acts.examples.hepmc3.Compression, args.hepmc3_compression
    )
    vtxGen = acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(
            args.beamspot_mean_mm[0] * u.mm,
            args.beamspot_mean_mm[1] * u.mm,
            args.beamspot_mean_mm[2] * u.mm,
            args.beamspot_mean_time_ns * u.ns,
        ),
        stddev=acts.Vector4(
            args.beamspot_sigma_mm[0] * u.mm,
            args.beamspot_sigma_mm[1] * u.mm,
            args.beamspot_sigma_mm[2] * u.mm,
            args.beamspot_sigma_time_ns * u.ns,
        ),
    )

    sequencer = acts.examples.Sequencer(
        events=args.events,
        numThreads=args.threads,
        logLevel=acts.logging.INFO,
    )
    runPythia8(
        outputDir=outputDir,
        outputRoot=not args.no_root,
        outputCsv=args.csv,
        hardProcess=hardProcess,
        outputHepMC3=args.hepmc3,
        hepmc3Compression=hepmc3Compression,
        vtxGen=vtxGen,
        s=sequencer,
    ).run()
