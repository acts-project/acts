"""Equivalence test for two Gaussian resolution-fit backends: ROOT's
`TH1::Fit` (via `ActsPlugins::RootHistogramFit`) and `acts.examples.scipy`'s
`curve_fit`-based one. A synthetic `IAlgorithm` writes tracks/particles/hits
straight to the whiteboard; the real `TrackTruthMatcher` then derives the
matching, and two `PythonTrackParameterPerformanceWriter`s differing only in
`fitFunction` see bit-identical histograms so only the fit itself can differ.
"""

import numpy as np
import pytest

import acts
import acts.examples
import acts.examples.scipy as acts_scipy
from acts.examples import PythonTrackParameterPerformanceWriter

u = acts.UnitConstants

pytestmark = pytest.mark.root


class _SyntheticTrackAlgorithm(acts.examples.IAlgorithm):
    """Writes `nTracks` synthetic tracks/particles/measurement-particle-map
    entries to the whiteboard every event: truth d0 = 0, fitted d0 = a
    residual drawn from `sampler` (a callable `rng -> float`).

    Deliberately does NOT write a TrackParticleMatching itself; instead it
    writes one source link per track plus a MeasurementParticlesMap entry, so
    the real TrackTruthMatcher (see `_run_backends`) derives the matching
    itself, as it would from real digitized hits.
    """

    def __init__(self, sampler, nTracks, seed):
        super().__init__(name="SyntheticTrackAlgorithm", level=acts.logging.WARNING)
        self._sampler = sampler
        self._nTracks = nTracks
        self._rng = np.random.default_rng(seed)

        self.outputTracks = acts.examples.WriteDataHandle(
            self, acts.examples.ConstTrackContainer, "OutputTracks"
        )
        self.outputTracks.initialize("tracks")
        self.outputParticles = acts.examples.WriteDataHandle(
            self, acts.examples.SimParticleContainer, "OutputParticles"
        )
        self.outputParticles.initialize("particles_selected")
        self.outputMeasurementParticlesMap = acts.examples.WriteDataHandle(
            self, acts.examples.MeasurementParticlesMap, "OutputMeasurementParticlesMap"
        )
        self.outputMeasurementParticlesMap.initialize("measurement_particles_map")

    def execute(self, context):
        tc = acts.examples.TrackContainer()
        particles = acts.examples.SimParticleContainer()
        measurementParticlesMap = acts.examples.MeasurementParticlesMap()

        surface = acts.Surface.createPerigee(acts.Vector3(0, 0, 0))
        # Identity covariance makes pull == residual (enough for this test).
        cov = acts.BoundMatrix.Identity()
        geoId = acts.GeometryIdentifier()

        for i in range(self._nTracks):
            barcode = acts.examples.SimBarcode()
            barcode.particle = i
            particle = acts.examples.SimParticle(barcode, acts.PdgParticle.eMuon)
            # Transverse direction, matching the track parameters below --
            # avoids a degenerate intersection with the perigee surface.
            particle.direction = acts.Vector3(1, 0, 0)
            particle.absoluteMomentum = 1.0 * u.GeV
            particles.insert(particle)

            # One measurement per track, own truth particle -> always Matched.
            hitIndex = i
            measurementParticlesMap.insert(hitIndex, barcode)

            residual = self._sampler(self._rng)
            track = tc.makeTrack()
            track.referenceSurface = surface
            track.parameters = acts.BoundVector(residual, 0.0, 0.0, np.pi / 2, 1.0, 0.0)
            track.covariance = cov
            track.particleHypothesis = acts.ParticleHypothesis.muon
            track.nMeasurements = 1

            state = track.appendTrackState()
            state.typeFlags.isMeasurement = True
            state.uncalibratedSourceLink = acts.examples.IndexSourceLink(
                geoId, hitIndex
            ).toSourceLink()

        self.outputTracks(context, tc.makeConst())
        self.outputParticles(context, particles)
        self.outputMeasurementParticlesMap(context, measurementParticlesMap)

        return acts.examples.ProcessCode.SUCCESS


def _small_res_plot_config():
    """Shrink Eta/Phi/Pt from their 40-bin defaults to 2 bins each -- every
    synthetic track lands in the same (eta, phi, pT) bin, so this just avoids
    fitting ~1600 empty slices per parameter for nothing.
    """
    cfg = acts.examples.ResPlotToolConfig()
    cfg.varBinning["Eta"] = acts.Axis.regular(2, -4.0, 4.0, "#eta")
    cfg.varBinning["Phi"] = acts.Axis.regular(2, -np.pi, np.pi, "#phi")
    cfg.varBinning["Pt"] = acts.Axis.regular(2, 0.0, 100.0, "pT [GeV/c]")
    return cfg


def _run_backends(sampler, nTracks, seed):
    """Run the synthetic algorithm + the real TrackTruthMatcher once, score
    the result with both fit backends, and return `{backend: histogram_dict}`.
    """
    import acts.examples.root as acts_root

    s = acts.examples.Sequencer(events=1, numThreads=1, logLevel=acts.logging.WARNING)
    s.addAlgorithm(_SyntheticTrackAlgorithm(sampler, nTracks, seed))
    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.WARNING,
            config=acts.examples.TrackTruthMatcher.Config(
                inputTracks="tracks",
                inputParticles="particles_selected",
                inputMeasurementParticlesMap="measurement_particles_map",
                outputTrackParticleMatching="track_particle_matching",
                outputParticleTrackMatching="particle_track_matching",
            ),
        )
    )

    writers = {}
    for backend, fitFn in [
        ("root", acts_root.makeRootHistogramFitFunction()),
        ("scipy", acts_scipy.makeScipyHistogramFitFunction()),
    ]:
        cfg = acts.examples.PythonTrackParameterPerformanceWriter.Config(
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            fitFunction=fitFn,
            resPlotToolConfig=_small_res_plot_config(),
        )
        writers[backend] = acts.examples.PythonTrackParameterPerformanceWriter(
            config=cfg, level=acts.logging.WARNING
        )
        s.addWriter(writers[backend])

    s.run()
    # histograms() re-runs fitProfiles() on every call; cache it once.
    return {backend: w.histograms() for backend, w in writers.items()}


def _fitted_bins(histograms, key, backend):
    """`(rootVals, otherVals, both)` for `key`, restricted to bins where both
    ROOT and `backend` succeeded (an unfitted bin has error == 0).
    """
    root = histograms["root"].get(key)
    other = histograms[backend].get(key)
    assert root is not None, f"ROOT produced no {key} (fit failed everywhere)"
    assert other is not None, f"{backend} produced no {key} (fit failed everywhere)"

    rootVals = np.asarray(root.histogram.values())
    rootErrs = np.asarray(root.histogram.errors())
    otherVals = np.asarray(other.histogram.values())
    otherErrs = np.asarray(other.histogram.errors())

    both = (rootErrs > 0) & (otherErrs > 0)
    assert (
        np.count_nonzero(both) >= 1
    ), f"no bin where both root and {backend} succeeded fitting {key}"
    return rootVals[both], otherVals[both], both


def _assert_backend_agrees(histograms, key, backend, rtol, atol):
    rootVals, otherVals, _ = _fitted_bins(histograms, key, backend)
    np.testing.assert_allclose(
        otherVals,
        rootVals,
        rtol=rtol,
        atol=atol,
        err_msg=f"{backend} vs root disagree on {key}",
    )


# Observed agreement is ~1e-5 to 1e-6 relative; 1e-3 leaves ample margin.
# resmean also gets an absolute floor since it's close to zero.
_RTOL = 1e-3
_MEAN_ATOL = 1e-3

_SCENARIOS = {
    "gaussian": lambda rng: rng.normal(0.0, 0.02),
    "gaussian_with_outliers": lambda rng: (
        rng.uniform(-0.5, 0.5) if rng.uniform() < 0.02 else rng.normal(0.0, 0.02)
    ),
}

# Fixed, not hash(scenario): Python randomizes string hashing per-process
# (PYTHONHASHSEED), so seeding off hash() would make this test's sample
# non-reproducible from run to run.
_SEEDS = {"gaussian": 2, "gaussian_with_outliers": 3}


@pytest.mark.parametrize("scenario", list(_SCENARIOS.keys()))
def test_fit_backends_agree(scenario):
    histograms = _run_backends(
        _SCENARIOS[scenario], nTracks=5000, seed=_SEEDS[scenario]
    )

    _assert_backend_agrees(
        histograms, "reswidth_d0_vs_eta", "scipy", rtol=_RTOL, atol=0.0
    )
    _assert_backend_agrees(
        histograms, "resmean_d0_vs_eta", "scipy", rtol=_RTOL, atol=_MEAN_ATOL
    )
