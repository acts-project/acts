"""
Examples:
  cosmic_particle_gun.py --events 1000 --generate-only
  cosmic_particle_gun.py --events 1000 --simulate --output cosmics
  cosmic_particle_gun.py --simulate --bz 2 --space-points measurements --ckf
  cosmic_particle_gun.py --events 10000 --simulate --bz 2 --seeding-output-only

"""

import argparse
import gc
import sys
import tempfile
from pathlib import Path

import numpy as np

# no __pycache__ next to the scripts (seeding_validation.py is imported from here)
sys.dont_write_bytecode = True

MUON_MASS = 0.1056584  # GeV

# space point class, surfaces of ODD:
# 17: pixel barrel, 16/18: pixel endcap, 24/29: strip barrel, 23/25/28/30: strip endcap
ODD_PIXEL_BARREL_VOLUMES = [17]
ODD_PIXEL_ENDCAP_VOLUMES = [16, 18]
ODD_STRIP_BARREL_VOLUMES = [24, 29]
ODD_STRIP_ENDCAP_VOLUMES = [23, 25, 28, 30]

# geometry from SpacePointMaker
ODD_2D_MEASUREMENT_VOLUMES = [16, 17, 18, 23, 24, 25]
ODD_LONG_STRIP_VOLUMES = [28, 29, 30]


def odd_possible_space_points(volumes, layers, sensitives):
    isLong = np.isin(volumes, ODD_LONG_STRIP_VOLUMES)
    sides = {}
    for v, l, s in zip(volumes[isLong], layers[isLong], sensitives[isLong]):
        sides.setdefault((int(v), int(l), (int(s) + 1) // 2), set()).add(int(s) % 2)
    return int((~isLong).sum()) + sum(len(x) == 2 for x in sides.values())


# cosmic muon functions: sample_energy, sample_cos_zenith, generate_batch, generate_accepted, make_muon_reader, inspect_input_particles
def sample_energy(rng, n, emin, emax, gamma):
    if abs(gamma - 1.0) < 1e-9:
        return emin * (emax / emin) ** rng.random(n)
    k = 1.0 - gamma
    return (emin**k + rng.random(n) * (emax**k - emin**k)) ** (1.0 / k)


def sample_cos_zenith(rng, n, ctcut, power):
    k = power + 2.0
    return (ctcut**k + rng.random(n) * (1.0 - ctcut**k)) ** (1.0 / k)


def generate_batch(rng, n, a):
    x0 = rng.uniform(a.xvert_low, a.xvert_high, n)
    z0 = rng.uniform(a.zvert_low, a.zvert_high, n)
    y0 = np.full(n, a.yvert_val)

    c = sample_cos_zenith(rng, n, a.ctcut, a.cos_power)
    s = np.sqrt(1.0 - c * c)
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    direction = np.stack([s * np.cos(az), -c, s * np.sin(az)], axis=1)

    energy = sample_energy(rng, n, a.emin, a.emax, a.gamma)
    p = np.sqrt(np.maximum(energy**2 - MUON_MASS**2, 0.0))

    plus = rng.random(n) < a.charge_ratio / (1.0 + a.charge_ratio)
    return {
        "pos": np.stack([x0, y0, z0], axis=1),
        "dir": direction,
        "p": p,
        "pdg": np.where(plus, -13, 13),  # mu+ = -13, mu- = +13
        "q": np.where(plus, 1.0, -1.0),
    }


# geometric filtering: entry_cylinder, entry_sphere, generate_accepted
def entry_cylinder(pos, d, radius, halfz):
    a = d[:, 0] ** 2 + d[:, 1] ** 2
    b = 2.0 * (pos[:, 0] * d[:, 0] + pos[:, 1] * d[:, 1])
    c = pos[:, 0] ** 2 + pos[:, 1] ** 2 - radius**2
    disc = b * b - 4.0 * a * c
    with np.errstate(invalid="ignore", divide="ignore"):
        t_side = (-b - np.sqrt(disc)) / (2.0 * a)
        z_side = pos[:, 2] + t_side * d[:, 2]
        ok_side = (disc >= 0) & (t_side >= 0) & (np.abs(z_side) <= halfz)
        best = np.where(ok_side, t_side, np.inf)
        for sign in (1.0, -1.0):  # caps at z = +-halfz
            t = (sign * halfz - pos[:, 2]) / d[:, 2]
            px = pos[:, 0] + t * d[:, 0]
            py = pos[:, 1] + t * d[:, 1]
            ok = (t >= 0) & (px * px + py * py <= radius**2)
            best = np.minimum(best, np.where(ok, t, np.inf))
    return np.isfinite(best), best


def entry_sphere(pos, d, ip, radius):
    w = np.asarray(ip)[None, :] - pos
    tc = np.einsum("ij,ij->i", w, d)
    dist2 = np.einsum("ij,ij->i", w, w) - tc * tc
    t = tc - np.sqrt(np.maximum(radius**2 - dist2, 0.0))
    return (dist2 <= radius**2) & (t >= 0), t


def generate_accepted(rng, a, n_needed, batch=200_000, min_acceptance=1e-5):
    chunks, got, tried = [], 0, 0
    while got < n_needed:
        if tried >= 50 * batch and got < min_acceptance * tried:
            raise RuntimeError(
                "Only %d of %d generated muons pass the filter (acceptance below "
                "%g): check the plane of generation and the volume of acceptance."
                % (got, tried, min_acceptance)
            )
        b = generate_batch(rng, batch, a)
        tried += batch
        if a.filter == "cylinder":
            ok, t = entry_cylinder(b["pos"], b["dir"], a.vol_radius, a.vol_halfz)
        elif a.filter == "sphere":
            ok, t = entry_sphere(b["pos"], b["dir"], a.ip, a.radius)
        else:
            ok, t = np.ones(batch, bool), np.zeros(batch)
        if not ok.any():
            continue
        # start slightly inside the volume so as not to land right on its boundary
        start = b["pos"][ok] + (t[ok] + a.entry_margin)[:, None] * b["dir"][ok]
        chunks.append(
            {
                "pos": start,
                "dir": b["dir"][ok],
                "p": b["p"][ok],
                "pdg": b["pdg"][ok],
                "q": b["q"][ok],
            }
        )
        got += int(ok.sum())
    muons = {k: np.concatenate([c[k] for c in chunks])[:n_needed] for k in chunks[0]}
    return muons, got / tried


def make_muon_reader(muons, per_event, outputParticles="particles_generated"):
    """Reader that gives the generated muons to the Sequencer from memory,
    per_event muons per event, instead of a CSV file per event.

    Each muon is its own primary vertex (vertexPrimary = 1, 2, ...,
    particle = 1), and its final state is its initial state, as
    CsvParticleReader does. The positions and momenta keep double precision
    (the CSV files had 4 decimals in mm and 6 in GeV).
    """
    import acts
    import acts.examples

    u = acts.UnitConstants

    class CosmicMuonReader(acts.examples.IReader):
        def __init__(self):
            acts.examples.IReader.__init__(self, "CosmicMuonReader", acts.logging.INFO)
            self.nEvents = len(muons["p"]) // per_event
            self.outputParticles = acts.examples.WriteDataHandle(
                self, acts.examples.SimParticleContainer, "OutputParticles"
            )
            self.outputParticles.initialize(outputParticles)

        def availableEvents(self):
            return (0, self.nEvents)

        def read(self, context):
            particles = acts.examples.SimParticleContainer()
            first = int(context.eventNumber) * per_event
            for j in range(per_event):
                i = first + j
                barcode = acts.examples.SimBarcode()
                barcode.vertexPrimary = j + 1
                barcode.particle = 1
                muon = acts.examples.SimParticle(
                    barcode,
                    (
                        acts.PdgParticle.eAntiMuon
                        if muons["pdg"][i] < 0
                        else acts.PdgParticle.eMuon
                    ),
                    float(muons["q"][i]) * u.e,
                    MUON_MASS * u.GeV,
                )
                x, y, z = muons["pos"][i]
                muon.fourPosition = acts.Vector4(x * u.mm, y * u.mm, z * u.mm, 0.0)
                muon.direction = acts.Vector3(*muons["dir"][i])
                muon.absoluteMomentum = float(muons["p"][i]) * u.GeV
                particles.insert(
                    acts.examples.SimParticle(muon.initialState, muon.initialState)
                )
            self.outputParticles(context, particles)
            return acts.examples.ProcessCode.SUCCESS

    return CosmicMuonReader()


def inspect_input_particles(path, max_events):
    """(events to run, particles per vertex)"""
    import uproot

    tree = uproot.open(path)["particles"]
    particle = tree["particle"].array(library="np")
    perVertex = max((int(p.max()) for p in particle if len(p) > 0), default=1)
    return min(max_events, tree.num_entries), perVertex


# cosmic seeding
class CosmicSeedingOutput:
    def __init__(self, path, rootFile, writers, effAlg):
        self.path = path
        self.rootFile = rootFile
        self.writers = writers
        self.effAlg = effAlg

    def finish(self):
        nOut = int(self.effAlg.nOutOfAcceptance)
        nOutMatched = int(self.effAlg.nOutOfAcceptanceMatched)
        for w in self.writers:
            w.flush()
        effWriter, seedWriter, eventWriter = self.writers
        self.close()

        import uproot

        tag = "[seeding]"
        f = uproot.open(self.path)
        eff = f["Efficiency"].arrays(
            ["event", "barcode", "matched", "nPossibleSP"], library="np"
        )
        events = f["Events"].arrays(["nSeeds", "nPureSeeds"], library="np")
        matched = eff["matched"].astype(bool)
        print(
            f"{tag} {effWriter.totalWritten} reconstructible particles, "
            f"{seedWriter.totalWritten} seeds, {eventWriter.totalWritten} "
            f"events -> {self.path}"
        )
        if len(matched) > 0:
            print(
                f"{tag} efficiency (pure seed per particle) = "
                f"{int(matched.sum())}/{len(matched)} = {100 * matched.mean():.1f}%"
            )
            seedable = eff["nPossibleSP"] >= 3
            if seedable.any() and not seedable.all():
                m = matched[seedable]
                print(
                    f"{tag} efficiency over particles with >= "
                    f"3 possible space points = "
                    f"{int(m.sum())}/{len(m)} = {100 * m.mean():.1f}%"
                )
        nSeeds = int(events["nSeeds"].sum())
        if nSeeds > 0:
            print(
                f"{tag} seeds per event = {nSeeds / len(events['nSeeds']):.1f}, "
                f"purity = {events['nPureSeeds'].sum() / nSeeds:.3f}"
            )
        if nOut > 0:
            print(
                f"{tag} {nOut} particles outside the angular acceptance "
                f"of the finder ({nOutMatched} with a pure seed anyway)"
            )
        if seedWriter.totalDropped > 0:
            print(
                f"{tag} WARNING: TTree 'Seeds' truncated by --max-seed-rows "
                f"({seedWriter.totalDropped} seeds not written)"
            )

    def close(self):
        self.effAlg = None
        if self.rootFile is not None:
            self.rootFile.close()
            self.rootFile = None


def add_cosmic_seeding(s, args, outputDir, trackingGeometry):
    import acts
    import acts.examples
    import uproot

    import seeding_validation as cse

    u = acts.UnitConstants

    effInputs = {}
    calibrate = (
        args.space_points == "measurements"
        and args.sp_selection == "all"
        and args.strip_formation == "direction"
    )
    seedSpacePoints = "spacepoints_calibrated" if calibrate else "spacepoints"
    if args.space_points == "simhits":
        s.addAlgorithm(
            acts.examples.SimHitSpacePointMaker(
                level=acts.logging.INFO,
                inputSimHits="simhits",
                outputSpacePoints="spacepoints",
                sigmaR=args.sigma_r * u.mm,
                sigmaZ=args.sigma_z * u.mm,
                # class per volume
                outputSpacePointClasses="spacepointclasses",
                pixelBarrelVolumes=ODD_PIXEL_BARREL_VOLUMES,
                pixelEndcapVolumes=ODD_PIXEL_ENDCAP_VOLUMES,
                stripBarrelVolumes=ODD_STRIP_BARREL_VOLUMES,
                stripEndcapVolumes=ODD_STRIP_ENDCAP_VOLUMES,
                # surface of each space point = surface of its SimHit
                outputSpacePointSurfaces="spacepointsurfaces",
            )
        )
    else:
        import acts.examples.json

        actsDir = Path(__file__).resolve().parents[3]
        if args.sp_selection == "acts":
            geometrySelection = acts.examples.json.readJsonGeometryList(
                str(actsDir / "Examples/Configs/odd-seeding-config.json")
            )
            stripGeometrySelection = []
            effInputs["reconstructibleGeometry"] = geometrySelection
        else:
            geometrySelection = [
                acts.GeometryIdentifier(volume=v) for v in ODD_2D_MEASUREMENT_VOLUMES
            ]
            stripGeometrySelection = acts.examples.json.readJsonGeometryList(
                str(actsDir / "Examples/Configs/odd-strip-spacepoint-selection.json")
            )
            # module of long strips (28, 29, 30) = 1,2,3; module of short strips (23, 24, 25) = 4,5,6
            effInputs["possibleSpacePoints"] = odd_possible_space_points
        s.addAlgorithm(
            acts.examples.SpacePointMaker(
                level=acts.logging.INFO,
                inputMeasurements="measurement_subset",
                outputSpacePoints="spacepoints",
                trackingGeometry=trackingGeometry,
                geometrySelection=geometrySelection,
                stripGeometrySelection=stripGeometrySelection,
                cosmicStripSpacePoints=args.strip_formation in ("cosmic", "direction"),
                cosmicStripKeepAllPairs=calibrate,
                stripVertex=acts.Vector3(*args.strip_vertex),
                stripLengthTolerance=args.strip_length_tolerance,
                outputSpacePointClasses="spacepointclasses",
                pixelBarrelVolumes=ODD_PIXEL_BARREL_VOLUMES,
                pixelEndcapVolumes=ODD_PIXEL_ENDCAP_VOLUMES,
                stripBarrelVolumes=ODD_STRIP_BARREL_VOLUMES,
                stripEndcapVolumes=ODD_STRIP_ENDCAP_VOLUMES,
                # surface of each space point = surface of its first measurement (strip) or of its SimHit (pixel)
                outputSpacePointSurfaces="spacepointsurfaces",
                inputMeasurementParticlesMap="measurement_particles_map",
                outputSpacePointParticlesMap="spacepoint_particles_map",
            )
        )
        effInputs.update(
            inputSpacePoints=seedSpacePoints,
            inputSpacePointParticlesMap="spacepoint_particles_map",
        )

    cuts = cse.make_cosmic_seed_cuts(args)
    # ODD disk: |z| < 3030 mm; ODD barrel: r < 1080 mm
    cuts.zMin = -args.z_max
    cuts.zMax = args.z_max
    cuts.useStripCalibration = calibrate
    cuts.stripCalibrationTolerance = args.strip_calibration_tolerance
    # surface veto, to avoid doublets on the same surface (strip or pixel) of the same module
    cuts.rejectSameSurface = not args.no_same_surface_veto
    s.addAlgorithm(
        acts.examples.CosmicTripletSeedingAlgorithm(
            level=acts.logging.INFO,
            inputSpacePoints="spacepoints",
            outputSeeds="seeds",
            cuts=cuts,
            inputSpacePointClasses="spacepointclasses",
            inputSpacePointSurfaces="spacepointsurfaces",
            outputSpacePoints="spacepoints_calibrated" if calibrate else "",
        )
    )

    if args.space_points == "measurements" and args.write_space_points:
        import acts.examples.root

        # space_points.root: SpacePointMaker positions;
        # space_points_calibrated.root: seed positions
        for name in dict.fromkeys(["spacepoints", seedSpacePoints]):
            s.addWriter(
                acts.examples.root.RootSpacePointWriter(
                    level=acts.logging.INFO,
                    inputSpacePoints=name,
                    inputSimHits="simhits",
                    inputMeasurementSimHitsMap="measurement_simhits_map",
                    inputMeasurementParticlesMap="measurement_particles_map",
                    trackingGeometry=trackingGeometry,
                    filePath=str(
                        outputDir
                        / (name.replace("spacepoints", "space_points") + ".root")
                    ),
                )
            )

    # each muon = its own primary vertex
    if args.input_particles is not None:
        perVertex = args.input_particles_per_vertex
    else:
        perVertex = 1 if args.mode == "cosmic" else args.gun_particles

    def add_efficiency(fileName, inputSeeds, name):
        path = outputDir / fileName
        rootFile = uproot.recreate(path)
        writers = cse.make_efficiency_writers(
            rootFile, args.flush_every, args.max_seed_rows or None
        )
        effAlg = cse.EfficiencyRootWriter(
            *writers,
            minHits=args.min_hits_reconstructible,
            minCosThetaY=cuts.minCosThetaY,
            keepOutOfAcceptance=args.keep_out_of_acceptance,
            barcodeOf=lambda pid: (pid.vertexPrimary - 1) * perVertex + pid.particle,
            inputSeeds=inputSeeds,
            name=name,
            level=acts.logging.INFO,
            **effInputs,
        )
        s.addAlgorithm(effAlg)
        return CosmicSeedingOutput(path, rootFile, writers, effAlg)

    return [add_efficiency("seeding_efficiency.root", "seeds", "EfficiencyRootWriter")]


def add_cosmic_ckf(s, args, outputDir, trackingGeometry, field):
    import acts
    import acts.examples
    from acts.examples.reconstruction import (
        AmbiguityResolutionConfig,
        CkfConfig,
        TrackSelectorConfig,
        addAmbiguityResolution,
        addCKFTracks,
    )
    from acts.examples.simulation import (
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )

    u = acts.UnitConstants

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-3.0, 3.0),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    s.addAlgorithm(
        acts.examples.TrackParamsEstimationAlgorithm(
            level=acts.logging.INFO,
            inputSeeds="seeds",
            outputTrackParameters="estimatedparameters",
            outputSeeds="estimatedseeds",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            initialSigmas=[
                1 * u.mm,
                1 * u.mm,
                1 * u.degree,
                1 * u.degree,
                0 * u.e / u.GeV,
                1 * u.ns,
            ],
            initialSigmaQoverPt=0.1 * u.e / u.GeV,
            initialSigmaPtRel=0.1,
            initialVarInflation=[1.0] * 6,
            particleHypothesis=acts.ParticleHypothesis.muon,
        )
    )

    trackSelector = None
    if args.ckf_track_selector == "odd":
        trackSelector = TrackSelectorConfig(
            pt=(0.0, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=7,
            maxHoles=2,
            maxOutliers=2,
        )
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        trackSelector,
        CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=2,
            seedDeduplication=True,
            stayOnSeed=True,
            pixelVolumes=[16, 17, 18],
            stripVolumes=[23, 24, 25],
            maxPixelHoles=1,
            maxStripHoles=2,
            constrainToVolumes=[
                2,  # beam pipe
                32,
                4,  # beam pip gap
                16,
                17,
                18,  # pixel
                20,  # PST
                23,
                24,
                25,  # short strip
                26,
                8,  # long strip gap
                28,
                29,
                30,  # long strip
            ],
            stopAtPerigee=False,
        ),
        outputDirRoot=outputDir,
        writeCovMat=True,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
        ),
        outputDirRoot=outputDir,
        writeCovMat=True,
    )


# ACTS chain: particle generation, simulation, seeding, CKF
def run_acts(args, muons, n_events):
    import acts
    import acts.examples
    import acts.examples.root
    from acts.examples import RandomNumbers, Sequencer
    from acts.examples.simulation import (
        EtaConfig,
        MomentumConfig,
        ParticleConfig,
        addDigitization,
        addFatras,
        addParticleGun,
    )

    u = acts.UnitConstants
    outputDir = Path(args.output)
    outputDir.mkdir(parents=True, exist_ok=True)
    rnd = RandomNumbers(seed=args.seed)

    s = Sequencer(events=n_events, numThreads=args.jobs, outputDir=str(outputDir))
    s.config.logLevel = acts.logging.INFO
    # with --seeding-output-only the ROOT files of the generation, simulation and digitization
    # (particles, particles_simulation, hits, measurements) go to a temporary directory, deleted at the end (cleanupTemporary)
    # the writers are not removed: the random numbers of an algorithm depend on its position in the sequence
    # (RandomNumbers::generateSeed), so without them the same --seed would simulate other events
    simOutputDir = outputDir
    if args.seeding_output_only:
        s.cosmicTemporaryDir = tempfile.TemporaryDirectory(prefix="cosmic_gun_")
        simOutputDir = Path(s.cosmicTemporaryDir.name)

    if args.input_particles is not None:
        # generated by RootParticleWriter
        s.addReader(
            acts.examples.root.RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(args.input_particles),
                outputParticles="particles_generated",
            )
        )
        s.addWhiteboardAlias("particles", "particles_generated")
        s.addWhiteboardAlias("particles_generated_selected", "particles_generated")
    elif args.mode == "cosmic":
        reader = make_muon_reader(muons, args.muons_per_event)
        s.addReader(reader)
        s.cosmicMuonReader = reader
        s.addWhiteboardAlias("particles", "particles_generated")
        s.addWhiteboardAlias("particles_generated_selected", "particles_generated")
        if args.print_particles:
            s.addAlgorithm(
                acts.examples.ParticlesPrinter(
                    level=acts.logging.INFO, inputParticles="particles_generated"
                )
            )
        s.addWriter(
            acts.examples.root.RootParticleWriter(
                level=acts.logging.INFO,
                inputParticles="particles_generated",
                filePath=str(simOutputDir / "particles.root"),
            )
        )
    else:  # particle_gun
        addParticleGun(
            s,
            outputDirRoot=simOutputDir,
            momentumConfig=MomentumConfig(
                args.gun_pt_range[0] * u.GeV,
                args.gun_pt_range[1] * u.GeV,
                transverse=True,
            ),
            etaConfig=EtaConfig(*args.gun_eta_range, uniform=True),
            particleConfig=ParticleConfig(
                args.gun_particles, acts.PdgParticle.eMuon, randomizeCharge=True
            ),
            multiplicity=args.gun_multiplicity,
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(50 * u.mm, 50 * u.mm, 50 * u.mm, 1 * u.ns),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
        )

    detector = None
    seeding = []  # CosmicSeedingOutput of the triplets
    if args.simulate:
        from acts.examples.odd import (
            getOpenDataDetector,
            getOpenDataDetectorDirectory,
        )

        odd_dir = args.odd_dir or getOpenDataDetectorDirectory()
        detector = getOpenDataDetector(odd_dir=odd_dir)
        trackingGeometry = detector.trackingGeometry()
        # field on/off
        field = acts.ConstantBField(acts.Vector3(0, 0, args.bz * u.T))

        actsDir = Path(__file__).resolve().parents[3]
        digi_config = args.digi_config or (
            actsDir / "Examples/Configs/odd-digi-smearing-config.json"
        )

        addFatras(s, trackingGeometry, field, rnd=rnd, outputDirRoot=simOutputDir)
        addDigitization(
            s,
            trackingGeometry,
            field,
            digiConfigFile=digi_config,
            outputDirRoot=simOutputDir,
            rnd=rnd,
        )
        if not args.no_seeding:
            seeding = add_cosmic_seeding(s, args, outputDir, trackingGeometry)
            if args.ckf:
                add_cosmic_ckf(s, args, outputDir, trackingGeometry, field)

    return s, detector, seeding


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--mode",
        choices=["cosmic", "particle_gun"],
        default="cosmic",
        help="cosmic = Athena scheme; particle_gun = standard ACTS gun",
    )
    p.add_argument(
        "--input-particles",
        type=Path,
        default=None,
        help="already generated particles.root",
    )
    p.add_argument("--events", type=int, default=10)
    p.add_argument("--jobs", type=int, default=-1, help="-1 = all threads")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", type=Path, default=Path.cwd())

    g = p.add_argument_group("CosmicGenerator (Athena property names)")
    g.add_argument("--emin", type=float, default=10.0, help="GeV")
    g.add_argument("--emax", type=float, default=1000.0, help="GeV")
    g.add_argument("--gamma", type=float, default=2.7, help="spectral index E^-gamma")
    g.add_argument(
        "--cos-power", type=float, default=2.0, help="vertical intensity ~ cos^n(theta)"
    )
    g.add_argument(
        "--ctcut",
        type=float,
        default=0.342,
        help="minimum cos(theta); 0.342 = 70 degrees from the vertical",
    )
    g.add_argument("--xvert-low", type=float, default=-6000.0, help="mm")
    g.add_argument("--xvert-high", type=float, default=6000.0, help="mm")
    g.add_argument("--zvert-low", type=float, default=-6000.0, help="mm")
    g.add_argument("--zvert-high", type=float, default=6000.0, help="mm")
    g.add_argument(
        "--yvert-val",
        type=float,
        default=1500.0,
        help="mm; height of the generation plane (must be > --vol-radius)",
    )
    g.add_argument(
        "--charge-ratio", type=float, default=1.27, help="mu+/mu-; 1.0 = 50/50"
    )

    f = p.add_argument_group("acceptance filters (pre-simulation)")
    f.add_argument(
        "--filter",
        choices=["cylinder", "sphere", "none"],
        default="cylinder",
        help="cylinder = crosses a volume; sphere = OptimizeForCavern",
    )
    f.add_argument(
        "--vol-radius",
        type=float,
        default=1050.0,
        help="mm; radius of the volume/entry (inside the geometry world)",
    )
    f.add_argument(
        "--vol-halfz", type=float, default=3000.0, help="mm; half-length in z"
    )
    f.add_argument(
        "--ip",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        metavar=("IPx", "IPy", "IPz"),
        help="mm; for --filter sphere",
    )
    f.add_argument(
        "--radius",
        type=float,
        default=1200.0,
        help="mm; radius of the sphere for --filter sphere",
    )
    f.add_argument(
        "--entry-margin",
        type=float,
        default=1.0,
        help="mm that the muon goes in past the volume boundary",
    )
    f.add_argument("--muons-per-event", type=int, default=1)
    f.add_argument(
        "--print-particles", action="store_true", help="prints the generated particles"
    )

    s = p.add_argument_group("simulation (ACTS)")
    s.add_argument(
        "--simulate",
        action="store_true",
        help="Fatras + digitization with the OpenDataDetector (ODD)",
    )
    s.add_argument(
        "--odd-dir",
        type=Path,
        default=None,
        help="ODD checkout; by default the ODD_PATH variable ",
    )
    s.add_argument(
        "--digi-config",
        type=Path,
        default=None,
        help="digitization config; by default "
        "Examples/Configs/odd-digi-smearing-config.json",
    )
    s.add_argument(
        "--bz", type=float, default=0.0, help="constant Bz [T]; 0 = field off"
    )
    s.add_argument(
        "--generate-only",
        action="store_true",
        help="only generates the muons and writes particles.root (one file, "
        "usable with --input-particles); no simulation",
    )
    s.add_argument(
        "--seeding-output-only",
        action="store_true",
        help="particles.root, particles_simulation.root, hits.root and "
        "measurements.root go to a temporary directory that is deleted at "
        "the end: with --simulate only seeding_efficiency.root (and the "
        "timing of the Sequencer) stays. Same events as without it",
    )

    q = p.add_argument_group("--mode particle_gun")
    q.add_argument("--gun-particles", type=int, default=1)
    q.add_argument("--gun-multiplicity", type=int, default=1)
    q.add_argument(
        "--gun-pt-range",
        nargs=2,
        type=float,
        default=[10.0, 100.0],
        metavar=("MIN", "MAX"),
        help="pT in GeV",
    )
    q.add_argument(
        "--gun-eta-range",
        nargs=2,
        type=float,
        default=[-3.0, 3.0],
        metavar=("MIN", "MAX"),
    )

    k = p.add_argument_group("cosmic seeding (with --simulate)")
    k.add_argument(
        "--no-seeding",
        action="store_true",
        help="only simulation and digitization, no seeding",
    )
    k.add_argument(
        "--space-points",
        choices=["simhits", "measurements"],
        default="simhits",
        help="simhits = SimHitSpacePointMaker (true position, "
        "as with the Athena hits); measurements = "
        "SpacePointMaker in cosmic mode on the digitization",
    )
    k.add_argument(
        "--sp-selection",
        choices=["acts", "all"],
        default="acts",
        help="with --space-points measurements: acts = like the "
        "ACTS seeding in the ODD (odd-seeding-config.json, "
        "pixels only, no strips; denominator in those same "
        "layers); all = pixels + short strips + pairs of "
        "long strips (--strip-formation)",
    )
    k.add_argument(
        "--strip-formation",
        choices=["direction", "cosmic", "vertex"],
        default="direction",
        help="with --sp-selection all, long strips: direction = "
        "all stereo pairs, and the seeder recomputes each "
        "one with the triplet direction (like the ITk strip "
        "seeding in Athena); cosmic = closest "
        "approach between the strips, no direction; vertex = "
        "constrained to the vertex --strip-vertex as in "
        "collisions",
    )
    k.add_argument(
        "--strip-calibration-tolerance",
        type=float,
        default=1.1,
        help="with --strip-formation direction: margin of the "
        "recomputed crossing on the strips, in half strips (1.1, "
        "the Athena one)",
    )
    k.add_argument(
        "--strip-vertex",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        metavar=("X", "Y", "Z"),
        help="mm; vertex of --strip-formation vertex. Very high "
        "in +y (e.g. 0 1e6 0) is equivalent to assuming vertical "
        "tracks",
    )
    k.add_argument(
        "--strip-length-tolerance",
        type=float,
        default=0.01,
        help="fraction of the strip length by which the crossing "
        "may fall outside it (SpacePointMaker)",
    )
    k.add_argument(
        "--write-space-points",
        action="store_true",
        help="with --space-points measurements: writes "
        "space_points.root (RootSpacePointWriter, residuals "
        "against the truth)",
    )
    k.add_argument(
        "--sigma-r",
        type=float,
        default=0.1,
        help="sigma of r [mm] (only --space-points simhits)",
    )
    k.add_argument(
        "--sigma-z",
        type=float,
        default=0.1,
        help="sigma of z [mm] (only --space-points simhits)",
    )
    k.add_argument("--max-kink", type=float, default=0.05, help="rad")
    k.add_argument("--max-residual", type=float, default=5.0, help="mm")
    k.add_argument(
        "--min-delta-y",
        type=float,
        default=5.0,
        help="mm; minimum |dy| per segment (Athena's mindRadius)",
    )
    k.add_argument(
        "--max-delta-y",
        type=float,
        default=2100.0,
        help="mm; maximum |dy| per segment (Athena's maxdRadius, "
        "1000 in the ID; default: ODD diameter)",
    )
    k.add_argument(
        "--min-cos-theta-y",
        type=float,
        default=None,
        help="minimum |dy|/|d|, in the finder and in the definition of "
        "reconstructible (CosmicSeedCuts default: 0.3)",
    )
    k.add_argument(
        "--z-max",
        type=float,
        default=3100.0,
        help="mm; maximum |z| of the three space points (zMin/zMax). "
        "Default: ODD envelope, whose long strip "
        "disks reach |z|=3025.5 mm",
    )
    # binning de los SP en el seeding
    # by default, the grid (CartesianSpacePointGrid) is defined only along the z-axis—using `--grid-bins 1 1 12`
    # and `--grid-neighbors 0 0 3`—which results in 12 bins, each with a size of (2 x --z-max) / 12 = 516.7 mm,
    # and 3 neighboring bins on each side
    k.add_argument(
        "--binning",
        choices=["athena", "grid"],
        default="grid",
        help="grid (default) = equidistant bins in x, y and z [mm] of the "
        "CartesianSpacePointGrid (--grid-bins, --grid-neighbors, "
        "--grid-x-range, --grid-y-range; z range = --z-max); "
        "athena = Athena's z bins (CosmicSeedCuts::zBinEdges / "
        "zBinNeighbors)",
    )
    # number of grid bins per axis (x, y, z)
    k.add_argument(
        "--grid-bins",
        type=int,
        nargs=3,
        default=[1, 1, 12],
        metavar=("NX", "NY", "NZ"),
        help="with --binning grid: number of bins in x, y and z",
    )
    # neighboring bins on each side per axis; 3 in z by default
    k.add_argument(
        "--grid-neighbors",
        type=int,
        nargs=3,
        default=[0, 0, 3],
        metavar=("NX", "NY", "NZ"),
        help="with --binning grid: neighbouring bins on each side per axis",
    )
    # x and y range of the grid [mm]
    k.add_argument(
        "--grid-x-range",
        type=float,
        nargs=2,
        default=[-1100.0, 1100.0],
        metavar=("MIN", "MAX"),
        help="mm; with --binning grid (ODD: r < 1080 mm)",
    )
    k.add_argument(
        "--grid-y-range",
        type=float,
        nargs=2,
        default=[-1100.0, 1100.0],
        metavar=("MIN", "MAX"),
        help="mm; with --binning grid",
    )
    k.add_argument(
        "--field-cuts",
        choices=["auto", "on", "off"],
        default="auto",
        help="cuts that take the curvature into account, the same "
        "as in run_athena_cosmics_chain.sh: --use-pt-cut "
        "--use-longitudinal-compatibility "
        "--use-pt-dependent-kink-residual "
        "--subdetector-angular-cuts. auto = active with field "
        "(|--bz| > 0): with the fixed kink, curved tracks "
        "lose seeds; without field there is no curvature",
    )
    k.add_argument("--use-pt-cut", action="store_true")
    k.add_argument("--use-longitudinal-compatibility", action="store_true")
    k.add_argument("--use-pt-dependent-kink-residual", action="store_true")
    k.add_argument(
        "--no-covariance-scale",
        action="store_true",
        help="do not apply Athena's errorsc={16,16,100,16}",
    )
    k.add_argument(
        "--no-same-surface-veto",
        action="store_true",
        help="do not discard the middle-bottom and middle-top links "
        "between two space points on the same surface "
        "(Athena's (*r)->sur() == sur0, "
        "CosmicSeedCuts::rejectSameSurface)",
    )
    k.add_argument(
        "--subdetector-angular-cuts",
        action="store_true",
        help="pixel: cos(dphi) > --pixel-min-cos-dphi instead of the kink",
    )
    k.add_argument("--pixel-min-cos-dphi", type=float, default=0.96)
    k.add_argument("--min-hits-reconstructible", type=int, default=3)
    k.add_argument("--keep-out-of-acceptance", action="store_true")
    k.add_argument(
        "--max-seed-rows",
        type=int,
        default=2000000,
        help="truncates the TTree 'Seeds' after these rows; 0 = no limit",
    )
    k.add_argument("--flush-every", type=int, default=20000)

    c = p.add_argument_group("track finding (CKF, like full_chain_odd.py)")
    c.add_argument(
        "--ckf",
        action="store_true",
        help="seeds -> parameter estimation -> CKF, with the "
        "full_chain_odd.py configuration. Requires "
        "--space-points measurements and |--bz| >= 0.1 T",
    )
    c.add_argument(
        "--ckf-track-selector",
        choices=["odd", "none"],
        default="none",
        help="none = no track selection (the whole "
        "cosmic is reconstructed); odd = TrackSelectorConfig from "
        "full_chain_odd.py, meant for collisions: "
        "|d0| < 4 mm with respect to (0,0,0) and <= 2 holes, which "
        "discard almost all cosmics",
    )
    return p


def main():
    args = build_parser().parse_args()
    if args.field_cuts == "on" or (args.field_cuts == "auto" and args.bz != 0.0):
        if args.bz == 0.0:
            raise SystemExit("--field-cuts on requires --bz != 0")
        args.use_pt_cut = True
        args.use_longitudinal_compatibility = True
        args.use_pt_dependent_kink_residual = True
        args.subdetector_angular_cuts = True
    if args.ckf:
        if not args.simulate or args.no_seeding or args.space_points != "measurements":
            raise SystemExit(
                "--ckf requires --simulate and --space-points "
                "measurements (without --no-seeding)"
            )
        if abs(args.bz) < 0.1:
            raise SystemExit(
                "--ckf requires |--bz| >= 0.1 T: the ACTS parameter "
                "estimation does not work without field"
            )
    if args.generate_only:
        if args.seeding_output_only:
            raise SystemExit(
                "--generate-only writes particles.root, "
                "--seeding-output-only does not"
            )
        args.simulate = False
    sequencer = None
    detector = None
    seeding = []
    try:
        muons, n_events = None, args.events
        if args.input_particles is not None:
            n_events, args.input_particles_per_vertex = inspect_input_particles(
                args.input_particles, args.events
            )
            print(
                f"[input] {n_events} events from {args.input_particles} "
                f"({args.input_particles_per_vertex} particles per vertex)"
            )
        elif args.mode == "cosmic":
            if args.filter == "cylinder" and args.yvert_val <= args.vol_radius:
                raise SystemExit("--yvert-val must be greater than --vol-radius")
            rng = np.random.default_rng(args.seed)
            n_muons = args.events * args.muons_per_event
            muons, eff = generate_accepted(rng, args, n_muons)
            print(
                f"[cosmic] acceptance efficiency = {100 * eff:.2f}% "
                f"-> {n_muons} accepted muons"
            )
        sequencer, detector, seeding = run_acts(args, muons, n_events)
        sequencer.run()
        for output in seeding:
            output.finish()
    finally:
        for output in seeding:
            output.close()
        temporary = getattr(sequencer, "cosmicTemporaryDir", None)
        del sequencer
        gc.collect()
        del detector
        gc.collect()
        if temporary is not None:
            temporary.cleanup()


if __name__ == "__main__":
    main()
