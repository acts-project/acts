"""Shared defaults, misalignment configuration, and CLI helpers (kfalignment workflow)."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

DEFAULT_NUM_EVENTS = 10000
DEFAULT_NUM_THREADS = -1
DEFAULT_INITIAL_VAR_INFLATION = [1] * 6

DEFAULT_SIMULATION_DIR = "simulation_output"
DEFAULT_ALIGNMENT_DIR = "alignment_output"
DEFAULT_RECONSTRUCTION_DIR = "reconstruction_output"


@dataclass
class MisalignmentConfig:
    """Misalignment target, scramble magnitudes, and alignment fit DoF mask."""

    target_volume: Union[int, list, dict] = -1
    target_layer: Union[int, list] = -1
    target_sensitive: Union[int, list] = -1
    target_extra: int = -1
    tx: int = 1
    ty: int = 1
    tz: int = 0
    rx: int = 0
    ry: int = 0
    rz: int = 1
    shift_mag_mm: float = 0.5
    rotation_mag_rad: float = 0.02
    seed: int = 42


DEFAULT_MISALIGNMENT_CONFIG = MisalignmentConfig()

_MISALIGN_TARGET_HELP = (
    "Which surfaces to misalign: '-1' or 'all' = all sensors; "
    "'VOL' = all layers in volume (e.g. 17); "
    "'VOL:LAY' = all sensors on that layer (e.g. 17:4); "
    "comma-separated for several (e.g. 17:8,24:4). "
    "Use --target=-1 (with '=') for all; no sensor-level selection."
)


@dataclass
class WorkflowPaths:
    """Resolved I/O paths for one workflow invocation."""

    srcdir: Path
    digi_config: Path
    simulation_dir: Path
    alignment_dir: Path
    reconstruction_dir: Path
    input_dir: Path
    output_dir: Path
    misalignment_dir: Path
    recon_mode: str = "nominal"


def _parse_one_misalignment_target(spec: str) -> Union[int, dict]:
    s = "".join(str(spec).split()).lower()
    if not s:
        raise ValueError("empty misalignment target token")

    if s in ("-1", "all", "*"):
        return -1

    if ":" in s:
        vol_s, lay_s = s.split(":", 1)
        if not vol_s or not lay_s or ":" in lay_s:
            raise ValueError(
                f"Invalid misalignment target token {spec!r}; "
                "expected '-1'/'all', 'VOL', or 'VOL:LAY'"
            )
        try:
            vol = int(vol_s)
            lay = int(lay_s)
        except ValueError as e:
            raise ValueError(
                f"Invalid misalignment target token {spec!r}; "
                "volume and layer must be integers"
            ) from e
        if vol == -1:
            raise ValueError(
                f"Invalid misalignment target token {spec!r}; "
                "use '-1' or 'all' alone for all sensors, not '-1:LAY'"
            )
        return {vol: {lay: -1}}

    try:
        vol = int(s)
    except ValueError as e:
        raise ValueError(
            f"Invalid misalignment target token {spec!r}; "
            "expected '-1'/'all', 'VOL', or 'VOL:LAY'"
        ) from e
    if vol == -1:
        return -1
    return {vol: -1}


def _merge_misalignment_targets(parts: list) -> Union[int, dict]:
    if any(p == -1 for p in parts):
        return -1

    merged: dict = {}
    for part in parts:
        if not isinstance(part, dict):
            raise ValueError(f"Unexpected misalignment target piece: {part!r}")
        for vol, layer_cfg in part.items():
            if vol not in merged:
                merged[vol] = layer_cfg
                continue
            cur = merged[vol]
            if cur == -1 or layer_cfg == -1:
                merged[vol] = -1
                continue
            if isinstance(cur, dict) and isinstance(layer_cfg, dict):
                cur.update(layer_cfg)
            elif isinstance(cur, dict):
                cur[layer_cfg] = -1
            elif isinstance(layer_cfg, dict):
                merged[vol] = {cur: -1, **layer_cfg}
            else:
                merged[vol] = {cur: -1, layer_cfg: -1}
    return merged


def parse_misalignment_target(spec: str) -> Union[int, dict]:
    raw = str(spec).strip()
    if not raw:
        raise ValueError("misalignment target is empty")

    tokens = [t for t in raw.split(",") if t.strip()]
    if not tokens:
        raise ValueError("misalignment target is empty")

    parts = [_parse_one_misalignment_target(t) for t in tokens]
    return _merge_misalignment_targets(parts)


def add_misalignment_dof_arguments(parser: argparse.ArgumentParser) -> None:
    d = DEFAULT_MISALIGNMENT_CONFIG
    parser.add_argument(
        "--tx",
        type=int,
        choices=[0, 1],
        default=d.tx,
        help="Enable local-x DoF (scramble + alignment fit, 0/1)",
    )
    parser.add_argument(
        "--ty",
        type=int,
        choices=[0, 1],
        default=d.ty,
        help="Enable local-y DoF (scramble + alignment fit, 0/1)",
    )
    parser.add_argument(
        "--tz",
        type=int,
        choices=[0, 1],
        default=d.tz,
        help="Enable local-z DoF for alignment fit only (never scrambled, 0/1)",
    )
    parser.add_argument(
        "--rx",
        type=int,
        choices=[0, 1],
        default=d.rx,
        help="Enable rx DoF for alignment fit only (never scrambled, 0/1)",
    )
    parser.add_argument(
        "--ry",
        type=int,
        choices=[0, 1],
        default=d.ry,
        help="Enable ry DoF for alignment fit only (never scrambled, 0/1)",
    )
    parser.add_argument(
        "--rz",
        type=int,
        choices=[0, 1],
        default=d.rz,
        help="Enable rz DoF (scramble + alignment fit, 0/1)",
    )
    parser.add_argument(
        "--shift-mag-mm",
        type=float,
        default=d.shift_mag_mm,
        help="Translation magnitude scale in mm",
    )
    parser.add_argument(
        "--rotation-mag-rad",
        type=float,
        default=d.rotation_mag_rad,
        help="Rotation magnitude scale in rad",
    )
    parser.add_argument(
        "--misalign-seed",
        type=int,
        default=d.seed,
        help=f"RNG seed for random misalignment scramble (default: {d.seed})",
    )


def misalignment_config_from_args(
    args, target: Optional[str] = None
) -> MisalignmentConfig:
    return misalignment_config_from_target(
        target if target is not None else getattr(args, "target", "-1"),
        tx=args.tx,
        ty=args.ty,
        tz=args.tz,
        rx=args.rx,
        ry=args.ry,
        rz=args.rz,
        shift_mag_mm=getattr(args, "shift_mag_mm", None),
        rotation_mag_rad=getattr(args, "rotation_mag_rad", None),
        seed=getattr(args, "misalign_seed", None),
    )


def misalignment_config_from_target(
    target: str = "-1",
    *,
    tx: Optional[int] = None,
    ty: Optional[int] = None,
    tz: Optional[int] = None,
    rx: Optional[int] = None,
    ry: Optional[int] = None,
    rz: Optional[int] = None,
    shift_mag_mm: Optional[float] = None,
    rotation_mag_rad: Optional[float] = None,
    seed: Optional[int] = None,
) -> MisalignmentConfig:
    d = DEFAULT_MISALIGNMENT_CONFIG
    return MisalignmentConfig(
        target_volume=parse_misalignment_target(target),
        target_layer=-1,
        target_sensitive=-1,
        target_extra=-1,
        tx=d.tx if tx is None else tx,
        ty=d.ty if ty is None else ty,
        tz=d.tz if tz is None else tz,
        rx=d.rx if rx is None else rx,
        ry=d.ry if ry is None else ry,
        rz=d.rz if rz is None else rz,
        shift_mag_mm=d.shift_mag_mm if shift_mag_mm is None else shift_mag_mm,
        rotation_mag_rad=(
            d.rotation_mag_rad if rotation_mag_rad is None else rotation_mag_rad
        ),
        seed=d.seed if seed is None else seed,
    )


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the top-level CLI parser (subcommands for each workflow stage)."""
    parser = argparse.ArgumentParser(
        description=(
            "ODD workflow: simulation, misalignment, alignment, and reconstruction. "
            "Run with no arguments for the full English usage tutorial."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python full_chain_odd_KFalignment.py                          # full tutorial
  python full_chain_odd_KFalignment.py full --target=17:8,24:4
  python full_chain_odd_KFalignment.py simulation
  python full_chain_odd_KFalignment.py misalignment --target=17:4
  python full_chain_odd_KFalignment.py alignment
  python full_chain_odd_KFalignment.py reconstruction --mode aligned
""",
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    full_parser = subparsers.add_parser(
        "full",
        help="Run complete workflow: simulation -> alignment -> reconstruction",
    )
    full_parser.add_argument(
        "--num-events",
        type=int,
        default=DEFAULT_NUM_EVENTS,
        help=f"Number of events (default: {DEFAULT_NUM_EVENTS})",
    )
    full_parser.add_argument(
        "--simulation-dir",
        type=Path,
        default=None,
        help=f"Simulation output directory (default: {DEFAULT_SIMULATION_DIR})",
    )
    full_parser.add_argument(
        "--alignment-dir",
        type=Path,
        default=None,
        help=f"Alignment output directory (default: {DEFAULT_ALIGNMENT_DIR})",
    )
    full_parser.add_argument(
        "--reconstruction-dir",
        type=Path,
        default=None,
        help=f"Reconstruction output directory (default: {DEFAULT_RECONSTRUCTION_DIR})",
    )
    full_parser.add_argument(
        "--target",
        type=str,
        default="-1",
        help=_MISALIGN_TARGET_HELP,
    )
    add_misalignment_dof_arguments(full_parser)

    sim_parser = subparsers.add_parser(
        "simulation",
        help="Run simulation on nominal ODD geometry",
    )
    sim_parser.add_argument(
        "--num-events",
        type=int,
        default=DEFAULT_NUM_EVENTS,
        help=f"Number of events (default: {DEFAULT_NUM_EVENTS})",
    )
    sim_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: {DEFAULT_SIMULATION_DIR})",
    )

    mis_parser = subparsers.add_parser(
        "misalignment",
        help="Apply random misalignment only (no alignment calibration)",
    )
    mis_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: {DEFAULT_ALIGNMENT_DIR})",
    )
    mis_parser.add_argument(
        "--target",
        type=str,
        default="-1",
        help=_MISALIGN_TARGET_HELP,
    )
    add_misalignment_dof_arguments(mis_parser)

    align_parser = subparsers.add_parser(
        "alignment",
        help="Run alignment calibration: read misalignment files, fit, write result",
    )
    align_parser.add_argument(
        "--num-events",
        type=int,
        default=DEFAULT_NUM_EVENTS,
        help=f"Number of events (default: {DEFAULT_NUM_EVENTS})",
    )
    align_parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help=f"Simulation input directory (default: {DEFAULT_SIMULATION_DIR})",
    )
    align_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Alignment output directory (default: {DEFAULT_ALIGNMENT_DIR})",
    )
    align_parser.add_argument(
        "--misalignment-dir",
        type=Path,
        default=None,
        help="Directory with misalignment files (default: same as --output-dir)",
    )
    add_misalignment_dof_arguments(align_parser)

    recon_parser = subparsers.add_parser(
        "reconstruction",
        help="Run reconstruction",
    )
    recon_parser.add_argument(
        "--mode",
        type=str,
        choices=["nominal", "misaligned", "aligned"],
        default="nominal",
        help="Reconstruction mode",
    )
    recon_parser.add_argument(
        "--num-events",
        type=int,
        default=DEFAULT_NUM_EVENTS,
        help=f"Number of events (default: {DEFAULT_NUM_EVENTS})",
    )
    recon_parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help=f"Input directory (default: {DEFAULT_SIMULATION_DIR})",
    )
    recon_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Reconstruction output base (default: {DEFAULT_RECONSTRUCTION_DIR}/<mode>)",
    )
    recon_parser.add_argument(
        "--alignment-dir",
        type=Path,
        default=None,
        help=f"Misalignment/aligned files (default: {DEFAULT_ALIGNMENT_DIR})",
    )

    return parser


def resolve_workflow_paths(
    args,
    command: str,
    srcdir: Path,
    script_dir: Path,
    digi_config: Path,
) -> WorkflowPaths:
    """Resolve all I/O paths from parsed args (centralized like full_chain_odd)."""
    cwd = Path.cwd()
    simulation_dir = (
        args.simulation_dir
        if command == "full" and getattr(args, "simulation_dir", None)
        else cwd / DEFAULT_SIMULATION_DIR
    )
    alignment_dir = (
        args.alignment_dir
        if command == "full" and getattr(args, "alignment_dir", None)
        else cwd / DEFAULT_ALIGNMENT_DIR
    )
    reconstruction_dir = (
        args.reconstruction_dir
        if command == "full" and getattr(args, "reconstruction_dir", None)
        else cwd / DEFAULT_RECONSTRUCTION_DIR
    )

    if command == "simulation":
        output_dir = args.output_dir or cwd / DEFAULT_SIMULATION_DIR
        return WorkflowPaths(
            srcdir=srcdir,
            digi_config=digi_config,
            simulation_dir=output_dir,
            alignment_dir=cwd / DEFAULT_ALIGNMENT_DIR,
            reconstruction_dir=cwd / DEFAULT_RECONSTRUCTION_DIR,
            input_dir=output_dir,
            output_dir=output_dir,
            misalignment_dir=cwd / DEFAULT_ALIGNMENT_DIR,
        )

    if command == "misalignment":
        output_dir = args.output_dir or cwd / DEFAULT_ALIGNMENT_DIR
        return WorkflowPaths(
            srcdir=srcdir,
            digi_config=digi_config,
            simulation_dir=cwd / DEFAULT_SIMULATION_DIR,
            alignment_dir=output_dir,
            reconstruction_dir=cwd / DEFAULT_RECONSTRUCTION_DIR,
            input_dir=cwd / DEFAULT_SIMULATION_DIR,
            output_dir=output_dir,
            misalignment_dir=output_dir,
        )

    if command == "alignment":
        input_dir = args.input_dir or cwd / DEFAULT_SIMULATION_DIR
        output_dir = args.output_dir or cwd / DEFAULT_ALIGNMENT_DIR
        misalignment_dir = args.misalignment_dir or output_dir
        return WorkflowPaths(
            srcdir=srcdir,
            digi_config=digi_config,
            simulation_dir=input_dir,
            alignment_dir=output_dir,
            reconstruction_dir=cwd / DEFAULT_RECONSTRUCTION_DIR,
            input_dir=input_dir,
            output_dir=output_dir,
            misalignment_dir=misalignment_dir,
        )

    if command == "reconstruction":
        input_dir = args.input_dir or cwd / DEFAULT_SIMULATION_DIR
        recon_base = args.output_dir or cwd / DEFAULT_RECONSTRUCTION_DIR
        mode = args.mode
        alignment_dir = args.alignment_dir or script_dir / DEFAULT_ALIGNMENT_DIR
        return WorkflowPaths(
            srcdir=srcdir,
            digi_config=digi_config,
            simulation_dir=input_dir,
            alignment_dir=alignment_dir,
            reconstruction_dir=recon_base,
            input_dir=input_dir,
            output_dir=recon_base / mode,
            misalignment_dir=alignment_dir,
            recon_mode=mode,
        )

    if command == "full":
        if args.simulation_dir:
            simulation_dir = args.simulation_dir
        if args.alignment_dir:
            alignment_dir = args.alignment_dir
        if args.reconstruction_dir:
            reconstruction_dir = args.reconstruction_dir
        return WorkflowPaths(
            srcdir=srcdir,
            digi_config=digi_config,
            simulation_dir=simulation_dir,
            alignment_dir=alignment_dir,
            reconstruction_dir=reconstruction_dir,
            input_dir=simulation_dir,
            output_dir=alignment_dir,
            misalignment_dir=alignment_dir,
        )

    raise ValueError(f"Unknown command: {command}")


def print_usage_tutorial() -> None:
    """Print a detailed English usage tutorial (no subcommand given)."""
    d = DEFAULT_MISALIGNMENT_CONFIG
    text = f"""
================================================================================
  full_chain_odd_KFalignment.py — ODD simulation / misalignment / alignment / reconstruction
================================================================================

USAGE
  python full_chain_odd_KFalignment.py <command> [options]

  Run without a command (this message) to see the full tutorial.
  For argparse help of one module:  python full_chain_odd_KFalignment.py <command> -h

--------------------------------------------------------------------------------
COMMANDS (modules)
--------------------------------------------------------------------------------

  full             End-to-end workflow:
                     1) simulation (nominal geometry)
                     2) misalignment (write misalignment files)
                     3) alignment fit (write aligned_result)
                     4) reconstruction in nominal / misaligned / aligned modes

  simulation       Generate particles, Fatras simulation, and digitization
                   on nominal ODD geometry.
                   Writes hits / measurements / particles under simulation_output.

  misalignment     Apply random local misalignment to selected sensitive surfaces.
                   Writes misalignment files (no track fit).

  alignment        Run track-based alignment calibration on misaligned geometry.
                   Requires prior simulation + misalignment outputs.
                   Writes aligned_result.txt (+ index map).

  reconstruction   Kalman reconstruction with a chosen geometry mode:
                     nominal | misaligned | aligned
                   Requires prior simulation (and misalignment/alignment files
                   for the non-nominal modes).

--------------------------------------------------------------------------------
SHARED DEFAULTS
--------------------------------------------------------------------------------

  Number of events (--num-events)     {DEFAULT_NUM_EVENTS}
  Sequencer threads                   {DEFAULT_NUM_THREADS}  (all CPU cores when -1)
  Alignment / scramble DoF mask       tx={d.tx} ty={d.ty} tz={d.tz}  rx={d.rx} ry={d.ry} rz={d.rz}
  Shift magnitude (--shift-mag-mm)    {d.shift_mag_mm} mm
  Rotation magnitude (--rotation-mag-rad)  {d.rotation_mag_rad} rad
  Misalignment RNG seed (--misalign-seed)  {d.seed}
  Default I/O directories:
    {DEFAULT_SIMULATION_DIR}/
    {DEFAULT_ALIGNMENT_DIR}/
    {DEFAULT_RECONSTRUCTION_DIR}/<mode>/

  Notes on DoFs:
    • Fit mask (--tx/--ty/--tz/--rx/--ry/--rz): which parameters alignment may fit.
    • Random scramble only ever injects tx, ty, rz (tz/rx/ry are never scrambled
      even if set to 1 for the fit).

--------------------------------------------------------------------------------
TYPICAL WORKFLOWS
--------------------------------------------------------------------------------

  A) One-shot end-to-end
     python full_chain_odd_KFalignment.py full --target=17:8,24:4 --num-events 10000

  B) Step by step
     python full_chain_odd_KFalignment.py simulation --num-events 10000
     python full_chain_odd_KFalignment.py misalignment --target=17:8,24:4
     python full_chain_odd_KFalignment.py alignment --num-events 10000
     python full_chain_odd_KFalignment.py reconstruction --mode nominal
     python full_chain_odd_KFalignment.py reconstruction --mode misaligned
     python full_chain_odd_KFalignment.py reconstruction --mode aligned

  C) Per-module help
     python full_chain_odd_KFalignment.py full -h
     python full_chain_odd_KFalignment.py misalignment -h
     python full_chain_odd_KFalignment.py alignment -h
     python full_chain_odd_KFalignment.py reconstruction -h

================================================================================
"""
    print(text.strip() + "\n")
