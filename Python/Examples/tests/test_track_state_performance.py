"""ROOT-free track-state residuals through the Python performance writer."""

import numpy as np
import pytest

import acts
import acts.examples as ex
from acts.examples.scipy import makeScipyHistogramFitFunction

pytestmark = pytest.mark.pypi


class SyntheticStates(ex.IAlgorithm):
    def __init__(self, truth, perigee=False):
        super().__init__(name="SyntheticStates", level=acts.logging.WARNING)
        self.perigee = perigee
        collections = {"tracks": ex.ConstTrackContainer}
        if truth:
            collections.update(
                particles=ex.SimParticleContainer,
                simhits=ex.SimHitContainer,
                measurement_particles=ex.MeasurementParticlesMap,
                measurement_simhits=ex.MeasurementSimHitsMap,
            )
        self.outputs = {}
        for name, kind in collections.items():
            handle = ex.WriteDataHandle(self, kind, name)
            handle.initialize(name)
            self.outputs[name] = handle

    def execute(self, ctx):
        tracks = ex.TrackContainer()
        particles = ex.SimParticleContainer()
        simhits = ex.SimHitContainer()
        measurement_particles = ex.MeasurementParticlesMap()
        measurement_simhits = ex.MeasurementSimHitsMap()
        barcode = ex.SimBarcode()
        barcode.particle = 1
        particle = ex.SimParticle(barcode, acts.PdgParticle.eMuon)
        particle.direction = acts.Vector3(1 / np.sqrt(2), 0, 1 / np.sqrt(2))
        particle.absoluteMomentum = np.sqrt(2)
        particles.insert(particle)
        track = tracks.makeTrack()
        track.particleHypothesis = acts.ParticleHypothesis.muon
        track.nMeasurements = 2
        if self.perigee:
            track.referenceSurface = acts.Surface.createPerigee(acts.Vector3(0, 0, 0))
            track.parameters = acts.BoundVector(
                0.125, 0, 0, np.pi / 4, -1 / np.sqrt(2), 0
            )
            track.covariance = acts.BoundMatrix.Identity()
        # State tests intentionally omit the track reference surface.
        for index, volume in enumerate((1, 2)):
            geo_id = acts.GeometryIdentifier(volume=volume, sensitive=1)
            surface = acts.Surface.createPlane(acts.Transform3.Identity(), None)
            surface.assignGeometryId(geo_id)
            state = track.appendTrackState()
            state.referenceSurface = surface
            state.typeFlags.setIsMeasurement()
            state.uncalibratedSourceLink = ex.IndexSourceLink(
                geo_id, index
            ).toSourceLink()
            state.smoothed = acts.BoundVector(
                0.125 + 0.1 * index, 0, 0, np.pi / 4, -1 / np.sqrt(2), 0
            )
            state.filtered = acts.BoundVector(
                -0.125, 0, 0, np.pi / 4, -1 / np.sqrt(2), 0
            )
            state.smoothedCovariance = acts.BoundMatrix.Identity()
            state.filteredCovariance = acts.BoundMatrix.Identity()
            state.allocateCalibrated(1)
            state.setProjectorSubspaceIndices([acts.eBoundLoc0])
            state.effectiveCalibrated = [0.025]
            state.effectiveCalibratedCovariance = [[2.0]]
            momentum = acts.Vector4(1, 0, 1, 2)
            simhits.insert(
                ex.SimHit(
                    geo_id, barcode, acts.Vector4(0, 0, 0, 0), momentum, momentum, index
                )
            )
            measurement_particles.insert(index, barcode)
            measurement_simhits.insert(index, index)
        values = dict(
            tracks=tracks.makeConst(),
            particles=particles,
            simhits=simhits,
            measurement_particles=measurement_particles,
            measurement_simhits=measurement_simhits,
        )
        for name, handle in self.outputs.items():
            handle(ctx, values[name])
        return ex.ProcessCode.SUCCESS


def writer_config(**kwargs):
    plots = ex.ResPlotToolConfig()
    for name in ("Eta", "Phi", "Pt"):
        plots.varBinning[name] = acts.Axis.regular(1, -4, 4, name)
    return ex.PythonTrackParameterPerformanceWriter.Config(
        inputTracks="tracks",
        resPlotToolConfig=plots,
        fitFunction=makeScipyHistogramFitFunction(),
        **kwargs,
    )


@pytest.mark.parametrize(
    "parameter_type, expected",
    [("Smoothed", 0.125), ("Filtered", -0.125), (None, 0.125)],
)
@pytest.mark.parametrize("select_volume", [False, True])
def test_track_state_truth_residuals(parameter_type, expected, select_volume):
    s = ex.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    s.addAlgorithm(SyntheticStates(truth=True))
    s.addAlgorithm(
        ex.TrackTruthMatcher(
            level=acts.logging.WARNING,
            inputTracks="tracks",
            inputParticles="particles",
            inputMeasurementParticlesMap="measurement_particles",
            outputTrackParticleMatching="matching",
            outputParticleTrackMatching="particle_matching",
        )
    )
    cfg = writer_config(
        inputParticles="particles",
        inputTrackParticleMatching="matching",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits",
        parameterSource=ex.TrackParameterSource.TrackState,
        parameterType=(
            getattr(ex.TrackParameterType, parameter_type) if parameter_type else None
        ),
        geometrySelection=[acts.GeometryIdentifier(volume=1)] if select_volume else [],
    )
    writer = ex.PythonTrackParameterPerformanceWriter(cfg, acts.logging.WARNING)
    s.addWriter(writer)
    s.run()
    histograms = writer.histograms()
    assert "res_d0" not in histograms
    hist = histograms["res_loc0"].histogram
    counts = np.asarray(hist.values())
    edges = np.asarray(hist.axis(0).edges)
    assert counts.sum() == (1 if select_volume else 2)
    assert counts[np.searchsorted(edges, expected, side="right") - 1] >= 1
    assert (
        np.asarray(histograms["resVsEta_loc0"].histogram.values()).sum() == counts.sum()
    )


def test_track_state_measurement_residuals_without_truth():
    s = ex.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    s.addAlgorithm(SyntheticStates(truth=False))
    cfg = writer_config(
        parameterSource=ex.TrackParameterSource.TrackState,
        reference=ex.TrackParameterReference.Measurement,
        parameterType=ex.TrackParameterType.Smoothed,
        geometrySelection=[acts.GeometryIdentifier(volume=1)],
    )
    writer = ex.PythonTrackParameterPerformanceWriter(cfg, acts.logging.WARNING)
    s.addWriter(writer)
    s.run()
    hist = writer.histograms()["res_loc0"].histogram
    counts = np.asarray(hist.values())
    edges = np.asarray(hist.axis(0).edges)
    assert counts.sum() == 1
    assert counts[np.searchsorted(edges, 0.125 - 0.025, side="right") - 1] == 1


@pytest.mark.parametrize(
    "missing, message",
    [
        ("inputSimHits", "Missing simulated hits input collection"),
        ("inputMeasurementSimHitsMap", "Missing measurement to simulated hits map"),
    ],
)
def test_track_state_truth_requires_hit_inputs(missing, message):
    cfg = writer_config(
        inputParticles="particles",
        inputTrackParticleMatching="matching",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits",
        parameterSource=ex.TrackParameterSource.TrackState,
    )
    setattr(cfg, missing, "")
    with pytest.raises(ValueError, match=message):
        ex.PythonTrackParameterPerformanceWriter(cfg, acts.logging.WARNING)


def test_perigee_defaults_need_no_hit_inputs():
    s = ex.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    s.addAlgorithm(SyntheticStates(truth=True, perigee=True))
    s.addAlgorithm(
        ex.TrackTruthMatcher(
            level=acts.logging.WARNING,
            inputTracks="tracks",
            inputParticles="particles",
            inputMeasurementParticlesMap="measurement_particles",
            outputTrackParticleMatching="matching",
            outputParticleTrackMatching="particle_matching",
        )
    )
    cfg = writer_config(
        inputParticles="particles", inputTrackParticleMatching="matching"
    )
    writer = ex.PythonTrackParameterPerformanceWriter(cfg, acts.logging.WARNING)
    s.addWriter(writer)
    s.run()
    histograms = writer.histograms()
    assert "res_loc0" not in histograms
    assert np.asarray(histograms["res_d0"].histogram.values()).sum() == 1
    assert writer.config.parameterSource == ex.TrackParameterSource.Track
    assert writer.config.reference == ex.TrackParameterReference.Truth
    assert writer.config.parameterType is None


def test_measurement_reference_requires_explicit_state_type():
    cfg = writer_config(
        parameterSource=ex.TrackParameterSource.TrackState,
        reference=ex.TrackParameterReference.Measurement,
    )
    with pytest.raises(ValueError, match="requires an explicit parameter type"):
        ex.PythonTrackParameterPerformanceWriter(cfg, acts.logging.WARNING)
