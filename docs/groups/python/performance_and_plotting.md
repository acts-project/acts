@defgroup python_performance_plotting Performance Evaluation and Plotting
@ingroup python_bindings
@brief Extracting efficiency/resolution numbers from a `Sequencer` run, and plotting them.

# Performance writers

ACTS has two output options for the same performance evaluation. Each pair below uses the same
collector and produces the same main histogram families. The Python writer returns histograms
through `.histograms()` after `s.run()`; the ROOT writer saves them in a ROOT file. ROOT output
also includes some summary objects and can include matching details.

| Evaluation | In-memory Python output (PyPI or source) | ROOT file output (ROOT-enabled source build) |
| :--- | :--- | :--- |
| Pattern recognition: efficiency, fake and duplicate tracks | `PythonPatternRecognitionPerformanceWriter` | `RootPatternRecognitionPerformanceWriter` |
| Track parameters: residuals, pulls, efficiency | `PythonTrackParameterPerformanceWriter` | `RootTrackParameterPerformanceWriter` |

The pairs take the same input collections and use similar configurations. For example, this
configures the Python pattern-recognition writer:

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

## Fit backend

The track-parameter writers fit residual and pull distributions to extract mean and width
profiles. The ROOT writer uses the ROOT fit backend. For the Python writer, set
`cfg.fitFunction = acts.examples.scipy.makeScipyHistogramFitFunction()` to use SciPy instead.

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
with matplotlib and mplhep:

```python
import matplotlib.pyplot as plt

histograms["trackeff_vs_pT"].plot()
plt.savefig("tracking_efficiency.svg")
```

The following plots illustrate the efficiency and residual views. They use example data, not
measured ACTS output; regenerate them with
`docs/examples/generate_python_performance_plots.py`.

![Illustrative tracking efficiency versus transverse momentum.](python/tracking_efficiency.svg){width=450px}
![Illustrative track-parameter residual distribution.](python/track_residual.svg){width=450px}

For a one-dimensional histogram, `boost_histogram.Histogram(histograms["res_d0"])` gives a
[boost-histogram](https://boost-histogram.readthedocs.io/) object for rebinning or other plotting
tools. With track-state evaluation, the first two parameter names become `loc0` and `loc1`
instead of `d0` and `z0`.

For geometry and track visualization (not histogram plotting), see
`acts.examples.visualization.PyVisualization2D` and `TrackVisualizerAlg`, and
`Examples/Scripts/generic_plotter.py` for the YAML-configured plotting tool used by physmon.

> [!note]
> `matplotlib`, `mplhep`, and `boost_histogram` are optional dependencies for this — install
> them alongside `pyacts` if you want to plot.
