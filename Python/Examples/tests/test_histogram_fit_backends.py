"""Equivalence test for two Gaussian resolution-fit backends: ROOT's
`TH1::Fit` (via `ActsPlugins::RootHistogramFit`) and a scipy `curve_fit`
callable. A synthetic `IAlgorithm` writes tracks/particles/hits straight to
the whiteboard; the real `TrackTruthMatcher` then derives the matching, and
two `PythonTrackParameterPerformanceWriter`s differing only in `fitFunction`
see bit-identical histograms so only the fit itself can differ.
"""

import numpy as np
import pytest

import acts
import acts.examples

try:
    from acts.examples import PythonTrackParameterPerformanceWriter
except ImportError:
    PythonTrackParameterPerformanceWriter = None

try:
    import acts.examples.root as acts_root
except ImportError:
    acts_root = None

u = acts.UnitConstants

pytestmark = [
    pytest.mark.root,
    pytest.mark.skipif(
        PythonTrackParameterPerformanceWriter is None or acts_root is None,
        reason="Python/ROOT performance writers not available",
    ),
]


def _gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def _gaussian_jac(x, amplitude, mean, sigma):
    """Analytic Jacobian of `_gaussian` w.r.t. (amplitude, mean, sigma).

    Routes curve_fit's underlying MINPACK call through `_lmder` (exact
    derivatives) instead of `_lmdif` (forward-difference approximation).
    """
    z = (x - mean) / sigma
    g = np.exp(-0.5 * z * z)
    return np.stack(
        [g, amplitude * g * z / sigma, amplitude * g * z * z / sigma], axis=-1
    )


def _scipy_gaussian_fit(hist, rng):
    """A ROOT-free Python fit backend using scipy.optimize.curve_fit.

    Matches `ActsExamples::HistogramFitFunction`'s signature. Drops empty
    bins rather than weighting them at sigma=1, mirroring ROOT's "SQ0". The
    analytic Jacobian plus a bounded maxfev/ftol/xtol keeps fitProfiles()
    fast versus MINPACK's defaults.
    """
    from scipy.optimize import curve_fit

    values = hist.histogram.values()
    edges = np.asarray(hist.histogram.axis(0).edges)
    centres = 0.5 * (edges[:-1] + edges[1:])

    if rng is not None:
        xMin, xMax = rng
        mask = (centres >= xMin) & (centres <= xMax)
        centres = centres[mask]
        values = values[mask]

    if np.count_nonzero(values) < 3 or values.sum() <= 0:
        return None

    mean0 = np.average(centres, weights=np.clip(values, 0, None))
    sigma0 = max(
        np.sqrt(np.average((centres - mean0) ** 2, weights=np.clip(values, 0, None))),
        1e-6,
    )
    amplitude0 = values.max()

    keep = values > 0
    fitCentres = centres[keep]
    fitValues = values[keep]
    errors = np.sqrt(fitValues)

    try:
        with np.errstate(all="ignore"):
            popt, pcov = curve_fit(
                _gaussian,
                fitCentres,
                fitValues,
                p0=[amplitude0, mean0, sigma0],
                sigma=errors,
                absolute_sigma=True,
                jac=_gaussian_jac,
                maxfev=200,
                ftol=1e-6,
                xtol=1e-6,
            )
    except RuntimeError:
        return None

    if not np.all(np.isfinite(pcov)):
        return None

    meanError = float(np.sqrt(pcov[1, 1]))
    sigmaError = float(np.sqrt(pcov[2, 2]))
    return (float(popt[1]), abs(float(popt[2])), meanError, sigmaError)


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
    cfg = acts_root.ResPlotToolConfig()
    cfg.varBinning["Eta"] = acts.Axis.regular(2, -4.0, 4.0, "#eta")
    cfg.varBinning["Phi"] = acts.Axis.regular(2, -np.pi, np.pi, "#phi")
    cfg.varBinning["Pt"] = acts.Axis.regular(2, 0.0, 100.0, "pT [GeV/c]")
    return cfg


def _run_backends(sampler, nTracks, seed):
    """Run the synthetic algorithm + the real TrackTruthMatcher once, score
    the result with both fit backends, and return `{backend: histogram_dict}`.
    """
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
        ("scipy", _scipy_gaussian_fit),
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

# "uniform" has no reswidth/resmean tolerance: a Gaussian fit to a flat-top
# distribution has no unique best fit, so backends legitimately disagree by
# up to several hundred percent. It's a pure existence/sanity check instead.
_SCENARIOS = {
    "uniform": lambda rng: rng.uniform(-0.05, 0.05),
    "gaussian": lambda rng: rng.normal(0.0, 0.02),
    "gaussian_with_outliers": lambda rng: (
        rng.uniform(-0.5, 0.5) if rng.uniform() < 0.02 else rng.normal(0.0, 0.02)
    ),
}

# Fixed, not hash(scenario): Python randomizes string hashing per-process
# (PYTHONHASHSEED), so seeding off hash() would make this test's sample
# non-reproducible from run to run.
_SEEDS = {"uniform": 1, "gaussian": 2, "gaussian_with_outliers": 3}


@pytest.mark.parametrize("scenario", list(_SCENARIOS.keys()))
def test_fit_backends_agree(scenario):
    pytest.importorskip("scipy")

    histograms = _run_backends(
        _SCENARIOS[scenario], nTracks=5000, seed=_SEEDS[scenario]
    )

    if scenario == "uniform":
        # Existence/sanity only (see the tolerance note above): a backend
        # is allowed to decline this ill-conditioned fit outright, but
        # whichever ones report a result must be well-formed.
        for backend in ["root", "scipy"]:
            hist = histograms[backend].get("reswidth_d0_vs_eta")
            assert hist is not None
            errs = np.asarray(hist.histogram.errors())
            vals = np.asarray(hist.histogram.values())
            fitted = errs > 0
            assert np.all(np.isfinite(vals[fitted]))
            assert np.all(vals[fitted] > 0)
        return

    _assert_backend_agrees(
        histograms, "reswidth_d0_vs_eta", "scipy", rtol=_RTOL, atol=0.0
    )
    _assert_backend_agrees(
        histograms, "resmean_d0_vs_eta", "scipy", rtol=_RTOL, atol=_MEAN_ATOL
    )
