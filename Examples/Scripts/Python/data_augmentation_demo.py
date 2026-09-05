from pathlib import Path
from typing import Optional

from pathlib import Path
from typing import Optional

import acts
import acts.examples
from acts.examples.root import (
    RootParticleReader,
    RootSimHitReader,
)
import uproot as ur
import awkward as ak
import numpy as np

u = acts.UnitConstants


def runTruthTracking(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputSimHitsPath: Optional[Path] = None,
    decorators=[],
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TrackSmearingSigmas,
        addTruthTrackingGsf,
    )
    from acts.examples.root import (
        RootParticleReader,
        RootTrackSummaryWriter,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)
    logger = acts.getDefaultLogger("GSF Example", acts.logging.INFO)

    if inputParticlePath is None:
        addParticleGun(
            s,
            ParticleConfig(num=1, pdg=acts.PdgParticle.eElectron, randomizeCharge=True),
            EtaConfig(-3.0, 3.0, uniform=True),
            MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0.015, 0.015, 55.0, 0),
            ),
            multiplicity=1,
            rnd=rnd,
            outputDirRoot=outputDir,
        )
    else:
        logger.info("Reading particles from {}", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )
        s.addWhiteboardAlias("particles", "particles_generated")

    if inputSimHitsPath is None:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
            enableInteractions=True,
        )
    else:
        logger.info("Reading hits from {}", inputSimHitsPath.resolve())
        s.addReader(
            RootSimHitReader(
                level=acts.logging.INFO,
                filePath=str(inputSimHitsPath.resolve()),
                outputSimHits="simhits",
            )
        )
        s.addWhiteboardAlias("particles_simulated_selected", "particles_generated")

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        trackSmearingSigmas=TrackSmearingSigmas(
            # zero everything so the GSF has a chance to find the measurements
            loc0=0,
            loc0PtA=0,
            loc0PtB=0,
            loc1=0,
            loc1PtA=0,
            loc1PtB=0,
            time=0,
            phi=0,
            theta=0,
            ptRel=0,
        ),
        particleHypothesis=acts.ParticleHypothesis.electron,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
    )

    addTruthTrackingGsf(
        s,
        trackingGeometry,
        field,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary.root"),
            writeGsfSpecific=True,
        )
    )

    return s


def readExampleRootData(outputDir: Path):
    outputDir = Path(outputDir)
    summary_tree = ur.open(outputDir / "tracksummary.root")["tracksummary"]
    particles_tree = ur.open(outputDir / "particles.root")["particles"]
    summary_fields = ["t_d0", "t_z0", "t_phi", "t_theta"]
    summary = summary_tree.arrays(summary_fields, library="ak")
    masks = []
    for field in summary_fields:
        nonempty = ak.to_numpy(ak.num(summary[field]) > 0)

        firsts = ak.firsts(summary[field])
        arr = ak.to_numpy(firsts)
        arr = np.asarray(arr).squeeze()
        # finite: numeric finite values (will be False for NaN/Inf)
        finite = np.isfinite(arr)

        combined_filter = np.logical_and(nonempty, finite)
        masks.append(combined_filter)

    particle_fields = ["vx", "vy", "vz", "px", "py", "pz", "q"]
    particles = particles_tree.arrays(particle_fields, library="ak")
    for field in particle_fields:
        nonempty = ak.to_numpy(ak.num(particles[field]) > 0)

        firsts = ak.firsts(particles[field])
        arr = ak.to_numpy(firsts)
        arr = np.asarray(arr).squeeze()
        finite = np.isfinite(arr)
        combined = np.logical_and(nonempty, finite)
        masks.append(combined)

    combined_mask = np.logical_and.reduce(masks)
    simulation_data = {
        "d0": ak.to_numpy(summary["t_d0"][combined_mask][:, 0]),
        "z0": ak.to_numpy(summary["t_z0"][combined_mask][:, 0]),
        "phi": ak.to_numpy(summary["t_phi"][combined_mask][:, 0]),
        "theta": ak.to_numpy(summary["t_theta"][combined_mask][:, 0]),
        "vx": ak.to_numpy(particles["vx"][combined_mask][:, 0]),
        "vy": ak.to_numpy(particles["vy"][combined_mask][:, 0]),
        "vz": ak.to_numpy(particles["vz"][combined_mask][:, 0]),
        "px": ak.to_numpy(particles["px"][combined_mask][:, 0]),
        "py": ak.to_numpy(particles["py"][combined_mask][:, 0]),
        "pz": ak.to_numpy(particles["pz"][combined_mask][:, 0]),
        "q": ak.to_numpy(particles["q"][combined_mask][:, 0]),
    }
    return simulation_data


def sampleUniformPointsIn3D(x_offset: float, y_offset: float, z_offset: float, n: int):
    if x_offset == 0 and y_offset == 0 and z_offset == 0:
        return np.zeros((n, 3))
    semi_axes = np.array([x_offset, y_offset, z_offset])
    points = np.empty((0, 3))
    # Expected acceptance rate of a candidate is pi/6 (volume of ellipsoid over
    # bounding box), so oversample accordingly to reduce the number of
    # rejection loop iterations.
    acceptance_rate = np.pi / 6
    while len(points) < n:
        remaining = n - len(points)
        n_candidates = int(remaining / acceptance_rate) + 10
        candidates = np.random.uniform(
            low=-semi_axes, high=semi_axes, size=(n_candidates, 3)
        )
        mask = np.sum((candidates / semi_axes) ** 2, axis=-1) <= 1
        points = np.vstack([points, candidates[mask]])
    return points[:n]


def voidTrackPropagation(
    beamspots,
    particle_info,
    return_beamspots=False,
):
    geo_context = acts.GeometryContext.dangerouslyDefaultConstruct()
    mag_field_context = acts.MagneticFieldContext()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    propagator = acts.EigenVoidPropagator(
        acts.EigenStepper(field), acts.VoidNavigator()
    )
    propagator_options = acts.PropagatorPlainOptions(geo_context, mag_field_context)
    # The particle hypothesis only affects the mass, which does not enter the
    # vacuum propagation used here (no material interactions), so a fixed pion
    # hypothesis is used regardless of the true particle species.
    particle_hypothesis = acts.ParticleHypothesis(211, 0.13957 * u.GeV, 1.0)

    n = particle_info.shape[0]
    new_poca_paramters = np.empty((n, 5))
    for i, truth_params_set in enumerate(particle_info):
        beamspot = beamspots[i]
        vtx = truth_params_set[:3]
        mom = truth_params_set[3:6]
        q = truth_params_set[6]

        pos4 = acts.Vector4(vtx[0], vtx[1], vtx[2], 0.0)
        qOverP = q / np.linalg.norm(mom)
        start = acts.BoundTrackParameters.createCurvilinear(
            pos4, acts.Vector3(*mom), qOverP, None, particle_hypothesis
        )
        target = acts.Surface.createPerigee(
            acts.Vector3(beamspot[0], beamspot[1], beamspot[2])
        )
        result = propagator.propagateToSurface(start, target, propagator_options)
        new_poca_paramters[i] = np.array(result.parameters)[:5]
    if return_beamspots:
        return new_poca_paramters, beamspots
    return new_poca_paramters


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    outputDir = Path.cwd()
    # Run a regular particle simulation
    runTruthTracking(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=outputDir,
    ).run()

    simulation_data = readExampleRootData(outputDir)

    particle_info = np.column_stack(
        [
            simulation_data["vx"],
            simulation_data["vy"],
            simulation_data["vz"],
            simulation_data["px"],
            simulation_data["py"],
            simulation_data["pz"],
            simulation_data["q"],
        ]
    )

    # Create new reference beamspots in an ellipsoid within the offset limits
    x_offset = 0.1
    y_offset = 1.1
    z_offset = 10
    beamspots = sampleUniformPointsIn3D(
        x_offset, y_offset, z_offset, len(particle_info)
    )

    # Calculate the PCA parameters for the new reference beamspots
    new_pcas = voidTrackPropagation(
        beamspots,
        particle_info,
        return_beamspots=False,
    )
    new_pcas = {
        "d0": new_pcas[:, 0],
        "z0": new_pcas[:, 1],
        "phi": new_pcas[:, 2],
        "theta": new_pcas[:, 3],
        "qOverP": new_pcas[:, 4],
    }
