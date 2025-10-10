"""HepMC3 utilities and normalization tools."""

from pathlib import Path

from acts._adapter import _patch_config
from acts import ActsPythonBindings

if not hasattr(ActsPythonBindings._examples, "_hepmc3"):
    raise ImportError("ActsPythonBindings._examples._hepmc3 not found")

_patch_config(ActsPythonBindings._examples._hepmc3)

from acts.ActsPythonBindings._examples._hepmc3 import *


def normalize(
    input_files,
    output_path=None,
    output_dir=".",
    output_prefix="events",
    events_per_file=10000,
    max_events=0,
    format=None,
    compression=None,
    compression_level=6,
    verbose=False,
):
    """
    Normalize and optionally chunk HepMC3 files.

    This function reads one or more HepMC3 files, normalizes event numbers,
    and writes them to output files. It can write to a single output file
    or chunk events into multiple files.

    Parameters
    ----------
    input_files : list[str|Path]
        List of input HepMC3 files to process
    output_path : str|Path, optional
        Single output file path. If specified, all events are written to this file.
        Format and compression are auto-detected from the filename.
        Mutually exclusive with chunking parameters.
    output_dir : str|Path, default="."
        Output directory for multi-file mode (ignored if output_path is set)
    output_prefix : str, default="events"
        Output file prefix for multi-file mode (ignored if output_path is set)
    events_per_file : int, default=10000
        Number of events per output file in multi-file mode (ignored if output_path is set)
    max_events : int, default=0
        Maximum number of events to process (0 = all events)
    format : Format, optional
        Output format (Format.ascii or Format.root). Auto-detected from output_path if provided.
        Defaults to Format.ascii in multi-file mode.
    compression : Compression, optional
        Compression type (Compression.none, .zlib, .lzma, .bzip2, .zstd).
        Auto-detected from output_path if provided. Defaults to Compression.none in multi-file mode.
    compression_level : int, default=6
        Compression level (0-19, higher = more compression)
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    HepMC3Normalizer.Result
        Result object containing:
        - numEvents: Number of events processed
        - outputFiles: List of output file paths created
        - totalInputSize: Total input size in bytes
        - totalOutputSize: Total output size in bytes
        - totalReadTime: Time spent reading (seconds)
        - totalWriteTime: Time spent writing (seconds)

    Examples
    --------
    >>> # Single output file with auto-detected compression
    >>> result = normalize(
    ...     input_files=["input1.hepmc3.gz", "input2.hepmc3"],
    ...     output_path="combined.hepmc3.zst",
    ...     verbose=True
    ... )

    >>> # Multi-file chunking
    >>> result = normalize(
    ...     input_files=["input.hepmc3"],
    ...     output_dir="normalized",
    ...     output_prefix="events",
    ...     events_per_file=1000,
    ...     compression=Compression.zstd,
    ...     compression_level=9
    ... )
    """
    config = HepMC3Normalizer.Config()

    # Convert input files to Path objects
    config.inputFiles = [Path(f) for f in input_files]

    # Single output mode
    if output_path is not None:
        output_path = Path(output_path)
        config.singleOutputPath = output_path
        # Auto-detect format and compression
        config.compression = (
            _detect_compression(output_path) if compression is None else compression
        )
        config.format = (
            formatFromFilename(str(output_path)) if format is None else format
        )
    else:
        # Multi-file mode
        config.outputDir = Path(output_dir)
        config.outputPrefix = output_prefix
        config.eventsPerFile = events_per_file
        config.compression = Compression.none if compression is None else compression
        config.format = Format.ascii if format is None else format

    config.maxEvents = max_events
    config.compressionLevel = compression_level
    config.verbose = verbose

    normalizer = HepMC3Normalizer(config)
    return normalizer.normalize()


def _detect_compression(path):
    """Detect compression type from file extension."""
    ext = path.suffix
    if ext == ".gz":
        return Compression.zlib
    elif ext == ".xz":
        return Compression.lzma
    elif ext == ".bz2":
        return Compression.bzip2
    elif ext == ".zst":
        return Compression.zstd
    return Compression.none
