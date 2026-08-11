"""Equivalence test for the three Gaussian resolution-fit backends: ROOT's
`TH1::Fit`, Core's own `gaussianHistogramFit`, and a scipy `curve_fit`
callable. A Python `IAlgorithm` writes a fixed set of synthetic
tracks/particles/measurement-particle-map straight to the whiteboard -- no
detector, no digitization -- with the fitted d0 residual drawn from an
engineered distribution (uniform, pure Gaussian, Gaussian with outliers). The
real `TrackTruthMatcher` algorithm then computes the track-particle matching
from that input, exactly as it would in a full reconstruction chain, rather
than a test fabricating the matching decision itself. Three
`PythonTrackFitterPerformanceWriter`s attached to the same whiteboard keys,
differing only in `fitFunction`, then see bit-identical histograms and only
the fit itself can differ.
"""

import numpy as np
import pytest

import acts
import acts.examples

try:
    from acts.examples import PythonTrackFitterPerformanceWriter
except ImportError:
    PythonTrackFitterPerformanceWriter = None

try:
    import acts.examples.root as acts_root
except ImportError:
    acts_root = None

u = acts.UnitConstants

pytestmark = [
    pytest.mark.root,
    pytest.mark.skipif(
        PythonTrackFitterPerformanceWriter is None or acts_root is None,
        reason="Python/ROOT performance writers not available",
    ),
]


def _gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def _scipy_gaussian_fit(hist, rng):
    """A ROOT-free Python fit backend using scipy.optimize.curve_fit.

    Matches ActsExamples::HistogramFitFunction's signature. Drops empty bins
    rather than weighting them at sigma=1, mirroring ROOT's "SQ0" / Core's
    gaussianHistogramFit, both of which give zero-content bins zero error and
    drop them from the least-squares sum.
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
                maxfev=10000,
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
    residual drawn from `sampler` (a callable `rng -> float`). All other
    track parameters are fixed so every track lands in the same (eta, phi,
    pT) bin.

    Deliberately does NOT write a TrackParticleMatching itself -- that would
    let the test assert the very match/fake/duplicate decision the real
    TrackTruthMatcher algorithm is responsible for making. Instead this
    writes one measurement-like source link per track plus a
    MeasurementParticlesMap entry tying it to the track's truth particle, and
    the real TrackTruthMatcher (run as a normal sequencer algorithm, see
    `_run_backends`) derives the matching from that, the same way it would
    from real digitized hits.
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
        # BoundMatrix has no Python setter beyond Zero()/Identity() -- Identity
        # makes pull == residual exactly (ResPlotTool divides by sqrt(cov_ii)),
        # which is enough for the resmean/reswidth comparison this test cares
        # about.
        cov = acts.BoundMatrix.Identity()
        geoId = acts.GeometryIdentifier()

        for i in range(self._nTracks):
            barcode = acts.examples.SimBarcode()
            barcode.particle = i
            particle = acts.examples.SimParticle(barcode, acts.PdgParticle.eMuon)
            # Transverse direction (theta = pi/2, eta = 0), matching the track
            # parameters below. A particle travelling along the perigee
            # surface's own axis would give a degenerate (parallel)
            # intersection, and ResPlotTool could not compute a truth
            # perigee parameter from it.
            particle.direction = acts.Vector3(1, 0, 0)
            particle.absoluteMomentum = 1.0 * u.GeV
            particles.insert(particle)

            # One measurement index per track, exclusively attributed to that
            # track's own truth particle -- TrackTruthMatcher's majority-hit
            # logic (1 of 1 hit, both reco- and truth-side) then always
            # yields a clean Matched classification.
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
    the result with all three fit backends, and return
    `{backend: histogram_dict}`.
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
        ("cpp", acts.examples.gaussianHistogramFit),
        ("scipy", _scipy_gaussian_fit),
    ]:
        cfg = acts.examples.PythonTrackFitterPerformanceWriter.Config(
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            fitFunction=fitFn,
            resPlotToolConfig=_small_res_plot_config(),
        )
        writers[backend] = acts.examples.PythonTrackFitterPerformanceWriter(
            config=cfg, level=acts.logging.WARNING
        )
        s.addWriter(writers[backend])

    s.run()
    # histograms() re-runs fitProfiles() on every call; cache it once.
    return {backend: w.histograms() for backend, w in writers.items()}


def _fitted_bins(histograms, key, backend):
    """`(rootVals, otherVals, both)` for `key`, restricted with a boolean mask
    to bins where both ROOT and `backend` succeeded (an unfitted
    Histogram bin is default-constructed at error == 0, which a genuine
    fitted width never is).
    """
    root = histograms["root"].get(key)
    other = histograms[backend].get(key)
    assert root is not None, f"ROOT produced no {key} (fit failed everywhere)"
    assert other is not None, f"{backend} produced no {key} (fit failed everywhere)"

    rootVals = np.asarray(root.values())
    rootErrs = np.asarray(root.errors())
    otherVals = np.asarray(other.values())
    otherErrs = np.asarray(other.errors())

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


# reswidth tolerances are tight for the two Gaussian-shaped scenarios
# (observed agreement is ~1e-5 to ~1e-6 relative; 1e-3 leaves ample margin
# against run-to-run float noise without masking a real regression). resmean
# additionally gets a small absolute floor: it is a residual mean genuinely
# close to zero, so its relative diff is dominated by near-zero-denominator
# bins and is not a meaningful check on its own -- the same caveat recorded
# for real reconstruction data in PROGRESS.md.
_RTOL = 1e-3
_MEAN_ATOL = 1e-3

# "uniform" has no reswidth/resmean tolerance at all: fitting a Gaussian to a
# flat-top distribution has no unique best fit, so LM/MINUIT/curve_fit
# legitimately settle at different points on a much flatter chi-square
# surface. Confirmed empirically that the disagreement is highly sample- and
# seed-dependent -- from a few percent up to several hundred percent for the
# exact same generative distribution -- so no fixed numeric tolerance would
# be both meaningful and stable. This is the same class of divergence as the
# "variable_bins" scenario excluded entirely (not given a loose tolerance)
# from Tests/UnitTests/Examples/Framework/GaussianHistogramFitRootBaselineTests.cpp.
# "uniform" is therefore a pure existence/sanity check: every backend must
# still produce a finite, positive-sigma fit, just not one that has to agree
# with the others.
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
        for backend in ["root", "cpp", "scipy"]:
            hist = histograms[backend].get("reswidth_d0_vs_eta")
            assert hist is not None
            errs = np.asarray(hist.errors())
            vals = np.asarray(hist.values())
            fitted = errs > 0
            assert np.all(np.isfinite(vals[fitted]))
            assert np.all(vals[fitted] > 0)
        return

    for backend in ["cpp", "scipy"]:
        _assert_backend_agrees(
            histograms, "reswidth_d0_vs_eta", backend, rtol=_RTOL, atol=0.0
        )
        _assert_backend_agrees(
            histograms, "resmean_d0_vs_eta", backend, rtol=_RTOL, atol=_MEAN_ATOL
        )
