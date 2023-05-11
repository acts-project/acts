#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

srcdir = Path(__file__).resolve().parent.parent.parent.parent
outputDir = Path.cwd()

# detector, trackingGeometry, _ = getOpenDataDetector()
detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()
field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

s = runTruthTrackingKalman(
    trackingGeometry,
    field,
    digiConfigFile=srcdir
    / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
    # "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
    outputDir=outputDir,
)

gsfOptions = {
    "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
    "maxComponents": 4,
    "abortOnError": False,
    "disableAllMaterialHandling": False,
    "finalReductionMethod": acts.examples.FinalReductionMethod.maxWeight,
    "weightCutoff": 1.0e-4,
    "level": acts.logging.INFO,
}

s.addAlgorithm(
    acts.examples.RefittingAlgorithm(
        acts.logging.INFO,
        inputTracks="kfTracks",
        outputTracks="gsfTracks",
        fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
    )
)

s.addAlgorithm(
    acts.examples.TracksToTrajectories(
        level=acts.logging.INFO,
        inputTracks="gsfTracks",
        outputTrajectories="gsfTrajectories",
    )
)

s.addWriter(
    acts.examples.TrackFitterPerformanceWriter(
        level=acts.logging.INFO,
        inputTrajectories="trajectories",
        inputParticles="truth_seeds_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        filePath=str(outputDir / "performance_refitter.root"),
    )
)

s.run()
