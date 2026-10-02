@defgroup python_performance_plotting Performance Evaluation and Plotting
@ingroup python_bindings
@brief Extracting efficiency/resolution numbers from a `Sequencer` run, and plotting them.

# Performance writers

ACTS has two output options for the same performance evaluation. Each pair uses the same
collector and produces the same main histogram families. The Python writer returns histograms
through `.histograms()` after `s.run()`; the ROOT writer saves them in a ROOT file. ROOT output
also includes some summary objects and can include matching details.

- **Pattern recognition** (efficiency, fake and duplicate tracks):
  `PythonPatternRecognitionPerformanceWriter` for in-memory output, or
  `RootPatternRecognitionPerformanceWriter` for a ROOT file.
- **Track parameters** (residuals, pulls, efficiency):
  `PythonTrackParameterPerformanceWriter` for in-memory output, or
  `RootTrackParameterPerformanceWriter` for a ROOT file.

The Python writers are available from PyPI or a source build. The ROOT writers require a
ROOT-enabled source build.

The pattern-recognition pair takes the same input collections and uses similar configurations.
With a `Sequencer` and truth matcher already configured, add the Python writer like this:

```python
cfg = acts.examples.PythonPatternRecognitionPerformanceWriter.Config()
cfg.inputTracks = "fitted_tracks"
cfg.inputParticles = "particles"
cfg.inputTrackParticleMatching = "track_particle_matching"
cfg.inputParticleTrackMatching = "particle_track_matching"
cfg.inputParticleMeasurementsMap = "particle_measurements_map"
writer = acts.examples.PythonPatternRecognitionPerformanceWriter(cfg, acts.logging.INFO)
s.addWriter(writer)
s.run()
histograms = writer.histograms()
```

For ROOT output, use `acts.examples.root.RootPatternRecognitionPerformanceWriter`, set the same
input collection fields on its `Config`, and set `filePath`. See the
[PyPI finding and fitting demo](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/pypi_finding_fitting_demo.py)
for the Python writers and the
[truth-tracking Kalman example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_kalman.py)
for ROOT output.

## Resolution fit backend

The track-parameter writers fit residual and pull distributions to extract mean and width
profiles. The ROOT writer uses the ROOT fit backend. For the Python writer, set
`cfg.fitFunction = acts.examples.scipy.makeScipyHistogramFitFunction()` to use SciPy instead;
install `scipy` and `numpy` separately. The ROOT track-parameter writer also supports evaluation
at individual track states and against calibrated measurements. The Python writer currently
evaluates track reference parameters against truth particles.

`TrackTruthMatcher(doubleMatching=True)` is the standard way to produce the
`inputTrackParticleMatching`/`inputParticleTrackMatching` collections these writers need.

## Available histograms

The histogram families are shared across output formats. The keys below are from the Python
writers' `.histograms()` dictionaries; the ROOT writers save corresponding histogram objects.
Names depend on configuration, and `pT`, `eta`, and `phi` variants are often available.

| Writer | Family | Example keys | Result |
| :--- | :--- | :--- | :--- |
| Pattern recognition | Tracking efficiency | `trackeff_vs_pT`, `trackeff_vs_eta` | `Efficiency1` |
| Pattern recognition | Fake and duplicate tracks | `fakeRatio_vs_pT`, `duplicationRatio_vs_eta` | `Efficiency1` |
| Pattern recognition | Candidate counts | `nRecoTracks_vs_pT`, `nFakeTracks_vs_eta` | `Histogram2` |
| Both | Track summary | `nMeasurements_vs_eta`, `nHoles_vs_pT` | `ProfileHistogram1` |
| Track parameters | Residuals and pulls | `res_d0`, `pull_d0` | `Histogram1` |
| Track parameters | Residuals versus kinematics | `resVsEta_d0`, `pullVsPt_phi` | `Histogram2` |
| Track parameters | Fitted mean and width | `resmean_d0_vs_eta`, `reswidth_d0_vs_eta` | `ProfileHistogram1` |

To get the complete list for *your* configuration after `s.run()`, including extra dimensions
and range-specific histograms:

```python
for name, histogram in sorted(writer.histograms().items()):
    print(f"{name:40} {type(histogram).__name__}")
```

## Plotting

One-dimensional `Histogram1`, `ProfileHistogram1`, and `Efficiency1` objects support `.plot()`
with matplotlib and mplhep. For example, `histograms["trackeff_vs_pT"].plot()` plots tracking
efficiency. This tested snippet uses a sample histogram:

@snippet{trimleft} examples/test_performance_and_plotting.py Plotting an ACTS histogram

The following plots illustrate the efficiency and residual views. They use randomly sampled
example counts, not measured ACTS output. The efficiency bars show binomial standard errors;
the residual-count bars show square-root count uncertainties. Regenerate them with
`docs/examples/generate_python_performance_plots.py`.

![Illustrative tracking efficiency versus transverse momentum.](python/tracking_efficiency.svg){width=450px}
![Illustrative track-parameter residual distribution.](python/track_residual.svg){width=450px}

ACTS `Histogram1` and `ProfileHistogram1` objects can also be converted to
[boost-histogram](https://boost-histogram.readthedocs.io/) objects for rebinning or other
plotting tools. The converted histogram can also be serialized with Python's `pickle`, so you
can save it and load it in a later analysis. For example, after the pattern-recognition writer
has run:

```python
import pickle

import boost_histogram as bh

measurements = bh.Histogram(histograms["nMeasurements_vs_eta"])
with open("n_measurements.pkl", "wb") as output:
    pickle.dump(measurements, output)
```

For an `Efficiency1`, convert its `.accepted` and `.total` histograms separately.

For geometry and track visualization (not histogram plotting), see
`acts.examples.visualization.PyVisualization2D` and `TrackVisualizerAlg`, and
`Examples/Scripts/generic_plotter.py` for the YAML-configured plotting tool used by physmon.

> [!note]
> `matplotlib`, `mplhep`, and `boost_histogram` are optional dependencies for this — install
> them alongside `pyacts` if you want to plot.
