#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants


def runRefittingKf(
    trackingGeometry,
    field,
    outputDir,
    s: acts.examples.Sequencer = None,
):
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    s = runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        # "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        outputDir=outputDir,
        s=s,
    )

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            level=acts.logging.VERBOSE,
            inputTracks="kf_tracks",
            outputTracks="kf_refit_tracks",
            fit=acts.examples.makeKfFitterFunction(trackingGeometry, field),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="truth_seeds_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_kf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_refitter.root"),
        )
    )

    return s


if __name__ == "__main__":
    outputDir = Path.cwd()

    # detector, trackingGeometry, decorators = getOpenDataDetector()
    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runRefittingKf(trackingGeometry, field, outputDir).run()
