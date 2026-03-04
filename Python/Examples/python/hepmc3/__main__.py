#!/usr/bin/env python3
"""Command-line interface for HepMC3 file normalization.

This module can be invoked as:
    python -m acts.examples.hepmc3 [options]
"""

import argparse
import json
import sys
from pathlib import Path

from . import (
    Compression,
    Format,
    availableCompressionModes,
    availableFormats,
    normalizeFiles,
)


def main(prog: str):

    parser = argparse.ArgumentParser(
        prog=prog,
        description="HepMC3 File Normalizer - Normalize and chunk HepMC files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Multi-file mode with chunking
  %(prog)s -i input1.hepmc3.gz input2.hepmc3 -n 1000 -c zstd -l 9
  %(prog)s -i file.root -o normalized/ -p events -n 5000

  # Single output mode (format/compression auto-detected)
  %(prog)s -i input1.hepmc3.gz input2.hepmc3 -S combined.hepmc3.zst
  %(prog)s -i input.hepmc3.gz -S output/events.hepmc3.gz
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        metavar="FILE",
        help="Input HepMC files",
    )

    parser.add_argument(
        "-S",
        "--single-output",
        metavar="FILE",
        help="Write all events to a single output file. Format and compression are detected from filename.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        metavar="DIR",
        help="Output directory (ignored with --single-output) [default: .]",
    )

    parser.add_argument(
        "-p",
        "--output-prefix",
        default="events",
        metavar="PREFIX",
        help="Output file prefix [default: events]",
    )

    parser.add_argument(
        "-n",
        "--events-per-file",
        type=int,
        default=10000,
        metavar="N",
        help="Number of events per output file (ignored with --single-output) [default: 10000]",
    )

    parser.add_argument(
        "-m",
        "--max-events",
        type=int,
        default=0,
        metavar="N",
        help="Maximum number of events to read (0 = all events) [default: 0]",
    )

    parser.add_argument(
        "-c",
        "--compression",
        choices=["none", "zlib", "lzma", "bzip2", "zstd"],
        default="none",
        help="Compression type (ignored with --single-output) [default: none]",
    )

    parser.add_argument(
        "-l",
        "--compression-level",
        type=int,
        default=6,
        metavar="LEVEL",
        help="Compression level (0-19, higher = more compression) [default: 6]",
    )

    parser.add_argument(
        "-f",
        "--format",
        choices=["ascii", "root"],
        default="ascii",
        help="Output format (ignored with --single-output) [default: ascii]",
    )

    parser.add_argument(
        "-j",
        "--json",
        action="store_true",
        help="Write JSON output with list of created files to stdout",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    parser.add_argument(
        "--list-compressions",
        action="store_true",
        help="List available compression modes and exit",
    )

    args = parser.parse_args()

    # Handle --list-compressions
    if args.list_compressions:
        print("Available compression modes:")
        for comp in availableCompressionModes():
            print(f"  {comp}")
        print("\nAvailable formats:")
        for fmt in availableFormats():
            print(f"  {fmt}")
        return 0

    # Convert string compression to enum
    compression_map = {
        "none": Compression.none,
        "zlib": Compression.zlib,
        "lzma": Compression.lzma,
        "bzip2": Compression.bzip2,
        "zstd": Compression.zstd,
    }

    # Convert string format to enum
    format_map = {
        "ascii": Format.ascii,
        "root": Format.root,
    }

    try:
        # Convert inputs to Path objects
        input_files = [Path(f) for f in args.input]
        single_output = Path(args.single_output) if args.single_output else None

        # Run normalization
        result = normalizeFiles(
            inputFiles=input_files,
            singleOutputPath=single_output,
            outputDir=Path(args.output_dir),
            outputPrefix=args.output_prefix,
            eventsPerFile=args.events_per_file,
            maxEvents=args.max_events,
            format=format_map.get(args.format),
            compression=compression_map.get(args.compression),
            compressionLevel=args.compression_level,
            verbose=args.verbose,
        )

        # Write JSON output if requested
        if args.json:
            output = {
                "num_events": result.numEvents,
                "num_files": len(result.outputFiles),
                "files": [
                    {
                        "path": str(Path(f).absolute()),
                        "size": Path(f).stat().st_size if Path(f).exists() else None,
                    }
                    for f in result.outputFiles
                ],
            }
            print(json.dumps(output, indent=2))

        return 0

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(prog="acts.examples.hepmc3"))
