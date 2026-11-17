import argparse
import gc
import math
import threading
from pathlib import Path

import numpy as np
import uproot

import acts
import acts.examples
import acts.examples.root
from acts.examples import Sequencer

u = acts.UnitConstants

# Athena barcode threshold (SIM_BARCODE_THRESHOLD)
DEFAULT_PRIMARY_BARCODE_MAX = 200000


def etaPhiCosThetaYFromExtremeHits(positions):
    top = positions[np.argmax(positions[:, 1])]
    bottom = positions[np.argmin(positions[:, 1])]
    d = bottom - top
    p = np.linalg.norm(d)
    if p == 0.0:
        return 0.0, 0.0, 0.0
    cosThetaY = abs(d[1]) / p
    phi = math.atan2(d[1], d[0])
    if p <= abs(d[2]):
        return 0.0, phi, cosThetaY
    return 0.5 * np.log((p + d[2]) / (p - d[2])), phi, cosThetaY


def etaPhiFromExtremeHitsVectorized(positions):
    n = positions.shape[0]
    rows = np.arange(n)
    top = positions[rows, np.argmax(positions[:, :, 1], axis=1)]
    bottom = positions[rows, np.argmin(positions[:, :, 1], axis=1)]
    d = bottom - top
    p = np.linalg.norm(d, axis=1)
    phi = np.where(p == 0.0, 0.0, np.arctan2(d[:, 1], d[:, 0]))
    dz = d[:, 2]
    ok = p > np.abs(dz)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.where(ok, 0.5 * np.log((p + dz) / (p - dz)), 0.0)
    return eta, phi


class IncrementalTreeWriter:
    def __init__(self, rootFile, treeName, columns, dtypes, flushEvery, maxRows=None):
        self.rootFile = rootFile
        self.treeName = treeName
        self.columns = columns
        self.dtypes = dtypes
        self.flushEvery = flushEvery
        self.maxRows = maxRows
        self.rows = []
        self.chunks = {c: [] for c in columns}
        self.nBuffered = 0
        self.totalWritten = 0
        self.totalDropped = 0
        rootFile.mktree(treeName, dict(zip(columns, dtypes)))

    def _room(self, n):
        if self.maxRows is None:
            return n
        return max(0, min(n, self.maxRows - self.totalWritten - self.nBuffered))

    def _rowsToChunk(self):
        if not self.rows:
            return
        for i, (c, dt) in enumerate(zip(self.columns, self.dtypes)):
            self.chunks[c].append(np.array([r[i] for r in self.rows], dtype=dt))
        self.rows = []

    def add(self, row):
        if self._room(1) == 0:
            self.totalDropped += 1
            return
        self.rows.append(row)
        self.nBuffered += 1
        if self.nBuffered >= self.flushEvery:
            self.flush()

    def addChunk(self, arrays):
        n = len(arrays[self.columns[0]])
        keep = self._room(n)
        self.totalDropped += n - keep
        if keep == 0:
            return
        self._rowsToChunk()
        for c, dt in zip(self.columns, self.dtypes):
            self.chunks[c].append(np.asarray(arrays[c][:keep], dtype=dt))
        self.nBuffered += keep
        if self.nBuffered >= self.flushEvery:
            self.flush()

    def flush(self):
        if self.nBuffered == 0:
            return
        self._rowsToChunk()
        arrays = {c: np.concatenate(self.chunks[c]) for c in self.columns}
        self.rootFile[self.treeName].extend(arrays)
        self.totalWritten += self.nBuffered
        self.chunks = {c: [] for c in self.columns}
        self.nBuffered = 0


class SpacePointRootWriter(acts.examples.IAlgorithm):
    def __init__(
        self, treeWriter, name="SpacePointRootWriter", level=acts.logging.INFO
    ):
        acts.examples.IAlgorithm.__init__(self, name, level)
        self.treeWriter = treeWriter
        self.inSpacePoints = acts.examples.ReadDataHandle(
            self, acts.SpacePointContainer, "spacepoints"
        )
        self.inSpacePoints.initialize("spacepoints")
        self._lock = threading.Lock()

    def execute(self, context):
        spacePoints = self.inSpacePoints(context.eventStore)
        n = len(spacePoints)
        if n == 0:
            return acts.examples.ProcessCode.SUCCESS
        with self._lock:
            self.treeWriter.addChunk(
                {
                    "event": np.full(n, int(context.eventNumber), dtype=np.int64),
                    "x": spacePoints.x,
                    "y": spacePoints.y,
                    "z": spacePoints.z,
                    "r": spacePoints.r,
                    "t": spacePoints.time,
                    "varianceR": spacePoints.varianceR,
                    "varianceZ": spacePoints.varianceZ,
                }
            )
        return acts.examples.ProcessCode.SUCCESS


class SeedRootWriter(acts.examples.IAlgorithm):
    def __init__(self, treeWriter, name="SeedRootWriter", level=acts.logging.INFO):
        acts.examples.IAlgorithm.__init__(self, name, level)
        self.treeWriter = treeWriter
        self.inSpacePoints = acts.examples.ReadDataHandle(
            self, acts.SpacePointContainer, "spacepoints"
        )
        self.inSpacePoints.initialize("spacepoints")
        self.inSeeds = acts.examples.ReadDataHandle(self, acts.SeedContainer, "seeds")
        self.inSeeds.initialize("seeds")
        self._lock = threading.Lock()

    def execute(self, context):
        seeds = self.inSeeds(context.eventStore)
        if len(seeds) == 0:
            return acts.examples.ProcessCode.SUCCESS
        spacePoints = self.inSpacePoints(context.eventStore)
        idx, quality, vertexZ = seeds.tripletArrays()
        keep = (idx >= 0).all(axis=1)
        idx, quality, vertexZ = idx[keep], quality[keep], vertexZ[keep]
        n = len(idx)
        if n == 0:
            return acts.examples.ProcessCode.SUCCESS
        x, y, z = spacePoints.x, spacePoints.y, spacePoints.z
        bot, mid, top = idx[:, 0], idx[:, 1], idx[:, 2]
        with self._lock:
            self.treeWriter.addChunk(
                {
                    "event": np.full(n, int(context.eventNumber), dtype=np.int64),
                    "x_bot": x[bot],
                    "y_bot": y[bot],
                    "z_bot": z[bot],
                    "x_mid": x[mid],
                    "y_mid": y[mid],
                    "z_mid": z[mid],
                    "x_top": x[top],
                    "y_top": y[top],
                    "z_top": z[top],
                    "quality": quality,
                    "vertexZ": vertexZ,
                }
            )
        return acts.examples.ProcessCode.SUCCESS


class EfficiencyRootWriter(acts.examples.IAlgorithm):

    def __init__(
        self,
        effWriter,
        seedWriter,
        eventWriter,
        minHits=3,
        minCosThetaY=0.0,
        keepOutOfAcceptance=False,
        primaryBarcodeMax=DEFAULT_PRIMARY_BARCODE_MAX,
        barcodeOf=None,
        inputSpacePoints=None,
        inputSpacePointParticlesMap=None,
        reconstructibleGeometry=None,
        possibleSpacePoints=None,
        inputSeeds="seeds",
        seedSize=3,
        name="EfficiencyRootWriter",
        level=acts.logging.INFO,
    ):
        acts.examples.IAlgorithm.__init__(self, name, level)
        self.effWriter = effWriter
        self.seedWriter = seedWriter
        self.eventWriter = eventWriter
        self.minHits = minHits
        self.seedSize = seedSize
        self.barcodeOf = barcodeOf or (lambda pid: pid.particle)
        self.minCosThetaY = minCosThetaY
        self.keepOutOfAcceptance = keepOutOfAcceptance
        self.nOutOfAcceptance = 0
        self.nOutOfAcceptanceMatched = 0
        self.primaryBarcodeMax = primaryBarcodeMax
        self.inSimHits = acts.examples.ReadDataHandle(
            self, acts.examples.SimHitContainer, "simhits"
        )
        self.inSimHits.initialize("simhits")
        self.inSeeds = acts.examples.ReadDataHandle(
            self, acts.SeedContainer, "InputSeeds"
        )
        self.inSeeds.initialize(inputSeeds)
        if (inputSpacePoints is None) != (inputSpacePointParticlesMap is None):
            raise ValueError(
                "inputSpacePoints and inputSpacePointParticlesMap " "go together"
            )
        self.inSpacePoints = None
        self.inSpacePointParticles = None
        if inputSpacePoints is not None:
            self.inSpacePoints = acts.examples.ReadDataHandle(
                self, acts.SpacePointContainer, "InputSpacePoints"
            )
            self.inSpacePoints.initialize(inputSpacePoints)
            self.inSpacePointParticles = acts.examples.ReadDataHandle(
                self,
                acts.examples.MeasurementParticlesMap,
                "InputSpacePointParticlesMap",
            )
            self.inSpacePointParticles.initialize(inputSpacePointParticlesMap)
        self.reconstructibleGeometry = None
        if reconstructibleGeometry is not None:
            self.reconstructibleGeometry = {
                (int(g.volume), int(g.layer)) for g in reconstructibleGeometry
            }
        self.possibleSpacePoints = possibleSpacePoints
        self._lock = threading.Lock()

    def execute(self, context):
        with self._lock:
            return self._executeLocked(context)

    def _executeLocked(self, context):
        simHits = self.inSimHits(context.eventStore)
        seeds = self.inSeeds(context.eventStore)
        event = int(context.eventNumber)

        nHits = len(simHits)
        if nHits == 0:
            self.eventWriter.add((event, 0, 0, 0, 0, 0))
            return acts.examples.ProcessCode.SUCCESS
        
        barcodes = np.empty(nHits, dtype=np.int64)
        positions = np.empty((nHits, 3), dtype=np.float64)
        momenta = np.empty((nHits, 3), dtype=np.float64)
        inGeometry = np.ones(nHits, dtype=bool)
        needGeometry = (
            self.reconstructibleGeometry is not None
            or self.possibleSpacePoints is not None
        )

        geometry = np.zeros((nHits, 3), dtype=np.int64)
        for i, h in enumerate(simHits):
            barcodes[i] = self.barcodeOf(h.particleId)
            pos4 = h.fourPosition
            positions[i] = (pos4[0], pos4[1], pos4[2])
            mom4 = h.momentum4Before
            momenta[i] = (mom4[0], mom4[1], mom4[2])
            if needGeometry:
                geoId = h.geometryId
                geometry[i] = (
                    int(geoId.volume),
                    int(geoId.layer),
                    int(geoId.sensitive),
                )
            if self.reconstructibleGeometry is not None:
                volume, layer = int(geometry[i, 0]), int(geometry[i, 1])
                inGeometry[i] = (volume, 0) in self.reconstructibleGeometry or (
                    volume,
                    layer,
                ) in self.reconstructibleGeometry

        spBarcodes, spPositions, nSpacePoints = barcodes, positions, nHits
        if self.inSpacePoints is not None:
            spacePoints = self.inSpacePoints(context.eventStore)
            spParticles = self.inSpacePointParticles(context.eventStore)
            nSpacePoints = len(spacePoints)
            spBarcodes = np.zeros(nSpacePoints, dtype=np.int64)
            spPositions = np.empty((nSpacePoints, 3), dtype=np.float64)
            if nSpacePoints > 0:
                spPositions[:, 0] = spacePoints.x
                spPositions[:, 1] = spacePoints.y
                spPositions[:, 2] = spacePoints.z
            nParticlesOf = np.zeros(nSpacePoints, dtype=np.int64)
            if len(spParticles) > 0:
                for i, pid in spParticles:
                    nParticlesOf[i] += 1
                    spBarcodes[i] = self.barcodeOf(pid)
            spBarcodes[nParticlesOf != 1] = 0

        matchedBarcodes = set()
        nSeeds = 0
        nPureSeeds = 0
        if len(seeds) > 0:
            idx = seeds.indexArrays(self.seedSize)[0]
            idx = idx[(idx >= 0).all(axis=1)]
            nSeeds = len(idx)
        if nSeeds > 0:
            bc = spBarcodes[idx]
            pure = (bc == bc[:, :1]).all(axis=1) & (bc[:, 0] > 0)
            nPureSeeds = int(pure.sum())
            matchedBarcodes = set(np.unique(bc[pure, 0]).tolist())
            eta, phi = etaPhiFromExtremeHitsVectorized(spPositions[idx])
            self.seedWriter.addChunk(
                {
                    "event": np.full(nSeeds, event, dtype=np.int64),
                    "eta": eta,
                    "phi": phi,
                    "pure": pure.astype(np.int32),
                }
            )

        if self.reconstructibleGeometry is not None:
            barcodes = barcodes[inGeometry]
            positions = positions[inGeometry]
            momenta = momenta[inGeometry]
            geometry = geometry[inGeometry]

        nParticles = 0
        nPrimaries = 0
        for bc in np.unique(barcodes):
            if bc <= 0:
                continue
            mask = barcodes == bc
            n = int(mask.sum())
            if n < self.minHits:
                continue
            hitPositions = positions[mask]
            eta, phi, cosThetaY = etaPhiCosThetaYFromExtremeHits(hitPositions)
            pTop = momenta[mask][np.argmax(hitPositions[:, 1])]
            pt = math.hypot(pTop[0], pTop[1])
            if pt == 0.0:
                pt = float("nan")
            isPrimary = int(bc <= self.primaryBarcodeMax)
            matched = bc in matchedBarcodes
            nPossible = n
            if self.possibleSpacePoints is not None:
                g = geometry[mask]
                nPossible = int(self.possibleSpacePoints(g[:, 0], g[:, 1], g[:, 2]))
            if cosThetaY < self.minCosThetaY:
                self.nOutOfAcceptance += 1
                self.nOutOfAcceptanceMatched += int(matched)
                if not self.keepOutOfAcceptance:
                    continue
            nParticles += 1
            nPrimaries += isPrimary
            self.effWriter.add(
                (
                    event,
                    int(bc),
                    isPrimary,
                    n,
                    nPossible,
                    float(eta),
                    float(phi),
                    float(pt),
                    float(cosThetaY),
                    int(matched),
                )
            )

        self.eventWriter.add((event, nSpacePoints, nParticles, nPrimaries, nSeeds, nPureSeeds))

        return acts.examples.ProcessCode.SUCCESS


def make_cosmic_seed_cuts(args):
    cuts = acts.examples.CosmicSeedCuts()
    cuts.maxKinkAngle = args.max_kink
    if args.min_cos_theta_y is not None:
        cuts.minCosThetaY = args.min_cos_theta_y
    cuts.maxResidual = args.max_residual
    cuts.minDeltaY = args.min_delta_y
    cuts.maxDeltaY = args.max_delta_y
    cuts.magneticFieldBz = args.bz
    cuts.usePtCut = args.use_pt_cut
    cuts.useLongitudinalCompatibility = args.use_longitudinal_compatibility
    cuts.usePtDependentKinkResidual = args.use_pt_dependent_kink_residual
    cuts.useCovarianceScale = not args.no_covariance_scale
    cuts.covarianceScale = [16.0, 16.0, 100.0, 16.0]
    cuts.useSubdetectorAngularCuts = args.subdetector_angular_cuts
    cuts.pixelMinCosDeltaPhi = args.pixel_min_cos_dphi
    return cuts


def make_efficiency_writers(rootFile, flushEvery, maxSeedRows=None):

    effTreeWriter = IncrementalTreeWriter(
        rootFile,
        "Efficiency",
        [
            "event",
            "barcode",
            "isPrimary",
            "nHits",
            "nPossibleSP",
            "eta",
            "phi",
            "pt",
            "cosThetaY",
            "matched",
        ],
        [
            np.int64,
            np.int64,
            np.int32,
            np.int32,
            np.int32,
            np.float64,
            np.float64,
            np.float64,
            np.float64,
            np.int32,
        ],
        flushEvery,
    )
    seedTreeWriter = IncrementalTreeWriter(
        rootFile,
        "Seeds",
        ["event", "eta", "phi", "pure"],
        [np.int64, np.float64, np.float64, np.int32],
        flushEvery,
        maxRows=maxSeedRows,
    )
    eventTreeWriter = IncrementalTreeWriter(
        rootFile,
        "Events",
        ["event", "nSpacePoints", "nParticles", "nPrimaries", "nSeeds", "nPureSeeds"],
        [np.int64, np.int64, np.int32, np.int32, np.int64, np.int64],
        flushEvery,
    )
    return effTreeWriter, seedTreeWriter, eventTreeWriter


def build_sequencer(args, rootFile, spRootFile=None):
    level = acts.logging.DEBUG if args.debug else acts.logging.INFO
    s = Sequencer(numThreads=args.jobs)
    s.config.logLevel = level

    s.addReader(
        acts.examples.root.RootSimHitReader(
            level=acts.logging.INFO,
            filePath=str(args.input),
            treeName="hits",
            outputSimHits="simhits",
        )
    )

    s.addAlgorithm(
        acts.examples.SimHitSpacePointMaker(
            level=level,
            inputSimHits="simhits",
            outputSpacePoints="spacepoints",
            sigmaR=args.sigma_r * u.mm,
            sigmaZ=args.sigma_z * u.mm,
            outputSpacePointClasses="spacepointclasses",
            stripMinRadius=args.strip_min_radius * u.mm,
            pixelBarrelMaxAbsZ=args.pixel_barrel_max_abs_z * u.mm,
            stripBarrelMaxAbsZ=args.strip_barrel_max_abs_z * u.mm,
        )
    )

    cuts = make_cosmic_seed_cuts(args)
    cuts.rejectSameSurface = False
    s.addAlgorithm(
        acts.examples.CosmicTripletSeedingAlgorithm(
            level=level,
            inputSpacePoints="spacepoints",
            outputSeeds="seeds",
            cuts=cuts,
            inputSpacePointClasses="spacepointclasses",
        )
    )

    spWriters = ()
    if spRootFile is not None:
        spTreeWriter = IncrementalTreeWriter(
            spRootFile,
            "SpacePoints",
            ["event", "x", "y", "z", "r", "t", "varianceR", "varianceZ"],
            [np.int64] + [np.float64] * 7,
            args.flush_every,
        )
        seedPositionTreeWriter = IncrementalTreeWriter(
            spRootFile,
            "Seeds",
            [
                "event",
                "x_bot",
                "y_bot",
                "z_bot",
                "x_mid",
                "y_mid",
                "z_mid",
                "x_top",
                "y_top",
                "z_top",
                "quality",
                "vertexZ",
            ],
            [np.int64] + [np.float64] * 11,
            args.flush_every,
        )
        s.addAlgorithm(SpacePointRootWriter(spTreeWriter, level=level))
        s.addAlgorithm(SeedRootWriter(seedPositionTreeWriter, level=level))
        spWriters = (spTreeWriter, seedPositionTreeWriter)

    effTreeWriter, seedTreeWriter, eventTreeWriter = make_efficiency_writers(
        rootFile, args.flush_every, args.max_seed_rows
    )
    effAlg = EfficiencyRootWriter(
        effTreeWriter,
        seedTreeWriter,
        eventTreeWriter,
        minHits=args.min_hits_reconstructible,
        minCosThetaY=cuts.minCosThetaY,
        keepOutOfAcceptance=args.keep_out_of_acceptance,
        primaryBarcodeMax=args.primary_barcode_max,
        level=level,
    )
    s.addAlgorithm(effAlg)
    return s, (effTreeWriter, seedTreeWriter, eventTreeWriter), effAlg, spWriters


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("input", type=Path, help="ROOT file from athena_to_acts_hits.py")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="output ROOT file, TTrees 'Efficiency', 'Seeds' and 'Events'",
    )
    p.add_argument(
        "--spacepoints-output",
        type=Path,
        default=None,
        help="additional ROOT file with the 'SpacePoints' and 'Seeds' TTrees "
        "(positions of each seed) from the same pass",
    )
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--sigma-r", type=float, default=0.1, help="sigma of r [mm]")
    p.add_argument("--sigma-z", type=float, default=0.1, help="sigma of z [mm]")
    p.add_argument(
        "--bz",
        type=float,
        default=2.0,
        help="T; assumed Bz",
    )
    p.add_argument("--max-kink", type=float, default=0.05)
    p.add_argument("--max-residual", type=float, default=5.0)
    p.add_argument(
        "--min-delta-y",
        type=float,
        default=5.0,
        help="mm; minimum |dy| per segment",
    )
    p.add_argument(
        "--max-delta-y",
        type=float,
        default=1000.0,
        help="mm; maximum |dy| per segment",
    )
    p.add_argument(
        "--min-cos-theta-y",
        type=float,
        default=None,
        help="minimum alignment |dy|/|d| with the y axis; applied "
        "to the finder AND to the definition of reconstructible. ",
    )
    p.add_argument(
        "--keep-out-of-acceptance",
        action="store_true",
        help="do not exclude from the denominator the particles with "
        "cosThetaY < minCosThetaY; measures acceptance x finder ",
    )
    p.add_argument("--use-pt-cut", action="store_true")
    p.add_argument("--use-longitudinal-compatibility", action="store_true")
    p.add_argument("--use-pt-dependent-kink-residual", action="store_true")
    p.add_argument(
        "--min-hits-reconstructible",
        type=int,
        default=3,
        help="minimum hits of a true particle in the event "
        "to count it as 'reconstructible'",
    )
    p.add_argument(
        "--primary-barcode-max",
        type=int,
        default=DEFAULT_PRIMARY_BARCODE_MAX,
        help="maximum barcode of a generator particle "
        "(isPrimary=1); above it are Geant4 secondaries",
    )
    p.add_argument(
        "--max-seed-rows",
        type=int,
        default=2000000,
        help="truncates the 'Seeds' TTree after this many rows (the "
        "per-event counts in 'Events' remain exact); 0 = no limit",
    )
    p.add_argument(
        "--no-covariance-scale",
        action="store_true",
        help="do NOT apply Athena's errorsc={16,16,100,16}; use the "
        "container variances unscaled",
    )
    p.add_argument(
        "--subdetector-angular-cuts",
        action="store_true",
        help="angular cut by type of the middle SP: strip uses only "
        "--max-kink; pixel uses cos(dphi) > --pixel-min-cos-dphi "
        "(xy plane) instead of the kink",
    )
    p.add_argument(
        "--pixel-min-cos-dphi",
        type=float,
        default=0.96,
        help="Athena's dFcut for pixel with --subdetector-angular-cuts "
        "(default 0.96, dphi < ~16 degrees)",
    )

    p.add_argument(
        "--strip-min-radius",
        type=float,
        default=350.0,
        help="r [mm] from which a hit is strip",
    )
    p.add_argument(
        "--pixel-barrel-max-abs-z",
        type=float,
        default=400.0,
        help="|z| [mm] barrel/endcap boundary in pixel",
    )
    p.add_argument(
        "--strip-barrel-max-abs-z",
        type=float,
        default=1400.0,
        help="|z| [mm] barrel/endcap boundary in strip",
    )
    p.add_argument("--flush-every", type=int, default=20000)
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()
    if args.max_seed_rows <= 0:
        args.max_seed_rows = None

    args.output.parent.mkdir(parents=True, exist_ok=True)
    s = None
    rootFile = uproot.recreate(args.output)
    spRootFile = None
    try:
        if args.spacepoints_output is not None:
            args.spacepoints_output.parent.mkdir(parents=True, exist_ok=True)
            spRootFile = uproot.recreate(args.spacepoints_output)
        s, writers, effAlg, spWriters = build_sequencer(args, rootFile, spRootFile)
        s.run()
        effTreeWriter, seedTreeWriter, eventTreeWriter = writers

        nOutOfAcceptance = int(effAlg.nOutOfAcceptance)
        nOutOfAcceptanceMatched = int(effAlg.nOutOfAcceptanceMatched)
        minCosThetaY = float(effAlg.minCosThetaY)
        effAlg = None
        for w in writers + spWriters:
            w.flush()
        print(
            f"[OK] {effTreeWriter.totalWritten} reconstructible particles, "
            f"{seedTreeWriter.totalWritten} seeds, "
            f"{eventTreeWriter.totalWritten} events -> {args.output}"
        )
        if spWriters:
            spTreeWriter, seedPositionTreeWriter = spWriters
            print(
                f"[OK] {spTreeWriter.totalWritten} space points, "
                f"{seedPositionTreeWriter.totalWritten} seeds -> "
                f"{args.spacepoints_output}"
            )

        if nOutOfAcceptance > 0:
            verb = "still counted" if args.keep_out_of_acceptance else "excluded"
            print(
                f"[INFO] {nOutOfAcceptance} particles outside the angular acceptance "
                f"of the finder ({verb}): cosThetaY < {minCosThetaY}, "
                f"which implies an upper limit on "
                f"|eta| = {math.acosh(1.0 / minCosThetaY):.3f}"
            )
            if nOutOfAcceptanceMatched > 0:
                print(
                    f"[NOTICE] {nOutOfAcceptanceMatched} of these "
                    f"({100.0 * nOutOfAcceptanceMatched / nOutOfAcceptance:.1f}%) "
                    f"had a 'pure' seed: the cut applied to the two extreme hits "
                    f"is stricter than the one used by the finder for the triplet. If this "
                    f"fraction is large, the acceptance proxy is not suitable for "
                    f"this dataset."
                )
        if seedTreeWriter.totalDropped > 0:
            print(
                    f"[NOTICE] The 'Seeds' TTree was truncated at "
                    f"{seedTreeWriter.totalWritten} rows: {seedTreeWriter.totalDropped} "
                    f"seeds were not written (--max-seed-rows). Global purity and "
                    f"seeds per event are derived from 'Events' and remain unaffected; "
                    f"purity vs. eta/phi, however, is affected (it is biased toward the initial events)."
            )
    finally:
        del s
        gc.collect()
        rootFile.close()
        if spRootFile is not None:
            spRootFile.close()


if __name__ == "__main__":
    main()
