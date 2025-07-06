#!/usr/bin/env python3

import argparse
from pathlib import Path
import re

import acts
import acts.examples
import acts.examples.hepmc3


def float_with_unit(s: str) -> float:
    m = re.match(r"(\d+\.?\d*)(.*)", s)
    assert m is not None, f"Invalid length string: {s}"

    num, unit = m.groups()

    scale = getattr(acts.UnitConstants, unit, 1)

    return float(num) * scale


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hard-scatter", "--hs", type=Path, help="Hard scatter file", required=True
    )
    parser.add_argument(
        "--hs-input-events",
        type=int,
        default=None,
        help="If None, the HepMC3 will read the entire file to determine the number of available events.",
    )
    parser.add_argument(
        "--pileup", "--pu", type=Path, help="Pileup file", required=True
    )
    parser.add_argument("--output", "-o", type=Path, help="Output file", required=True)
    parser.add_argument("--force", "-f", action="store_true", help="Force overwrite")
    parser.add_argument(
        "--pileup-multiplicity",
        "--npu",
        type=int,
        help="Pileup multiplicity",
        default=1,
    )

    parser.add_argument(
        "--vtx-pos",
        "--vp",
        type=float_with_unit,
        help="Vertex position, use suffix like mm for units",
        default=[0.0, 0.0, 0.0, 0.0],
        nargs=4,
    )

    parser.add_argument(
        "--vtx-stddev",
        "--vs",
        type=float_with_unit,
        help="Vertex stddev, use suffix like mm for units",
        default=[0.0, 0.0, 0.0, 0.0],
        nargs=4,
    )

    parser.add_argument("--jobs", "-j", type=int, help="Number of jobs", default=-1)

    args = parser.parse_args()

    if not args.hard_scatter.exists():
        raise FileNotFoundError(f"Hard scatter file {args.hard_scatter} does not exist")
    if not args.pileup.exists():
        raise FileNotFoundError(f"Pileup file {args.pileup} does not exist")

    extensions = {
        acts.examples.hepmc3.compressionExtension(c): c
        for c in acts.examples.hepmc3.availableCompressionModes()
        if c != acts.examples.hepmc3.Compression.none
    }

    compression = extensions.get(
        args.output.suffix, acts.examples.hepmc3.Compression.none
    )

    if compression != acts.examples.hepmc3.Compression.none:
        output_path = Path(args.output.stem)
    else:
        output_path = args.output

    if output_path.exists():
        if args.force:
            output_path.unlink()
        else:
            raise FileExistsError(
                f"Output file {args.output} already exists, run with --force to overwrite"
            )

    s = acts.examples.Sequencer(numThreads=args.jobs, logLevel=acts.logging.WARNING)

    rng = acts.examples.RandomNumbers(seed=42)
    s.addReader(
        acts.examples.hepmc3.HepMC3Reader(
            inputPaths=[
                (args.hard_scatter, 1),
                (args.pileup, args.pileup_multiplicity),
            ],
            level=acts.logging.INFO,
            outputEvent="hepmc3_event",
            checkEventNumber=False,  # This is not generally guaranteed for arbitrary inputs
            numEvents=args.hs_input_events,
            randomNumbers=rng,
            vertexGenerator=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    args.vtx_stddev[0],
                    args.vtx_stddev[1],
                    args.vtx_stddev[2],
                    args.vtx_stddev[3],
                ),
                mean=acts.Vector4(
                    args.vtx_pos[0],
                    args.vtx_pos[1],
                    args.vtx_pos[2],
                    args.vtx_pos[3],
                ),
            ),
        )
    )

    s.addWriter(
        acts.examples.hepmc3.HepMC3Writer(
            inputEvent="hepmc3_event",
            outputPath=output_path,
            level=acts.logging.INFO,
            compression=compression,
            writeEventsInOrder=False,
        )
    )

    s.run()


if __name__ == "__main__":
    main()
