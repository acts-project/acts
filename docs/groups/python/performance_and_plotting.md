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

@snippet{trimleft} "../Examples/Scripts/Python/pypi_finding_fitting_demo.py" Python pattern-recognition performance writer

The [PyPI demo test](https://github.com/acts-project/acts/blob/main/Python/Examples/tests/test_examples.py)
executes this example. After `s.run()`, retrieve the results with
`histograms = perfWriterFinder.histograms()`.

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
install `scipy` and `numpy` separately. Both track-parameter writers support evaluation at
individual track states and against calibrated measurements.

`TrackTruthMatcher(doubleMatching=True)` is the standard way to produce the
`inputTrackParticleMatching`/`inputParticleTrackMatching` collections these writers need.

## Track-state residuals

By default, track reference parameters are compared with truth particles. To evaluate
smoothed parameters at individual measurement surfaces instead, configure the Python writer:

```python
cfg = acts.examples.PythonTrackParameterPerformanceWriter.Config(
    inputTracks="tracks",
    inputParticles="particles",
    inputTrackParticleMatching="track_particle_matching",
    parameterSource=acts.examples.TrackParameterSource.TrackState,
    parameterType=acts.examples.TrackParameterType.Smoothed,
    inputSimHits="simhits",
    inputMeasurementSimHitsMap="measurement_simhits_map",
    fitFunction=acts.examples.scipy.makeScipyHistogramFitFunction(),
)
writer = acts.examples.PythonTrackParameterPerformanceWriter(cfg, acts.logging.INFO)
s.addWriter(writer)
```

The residual is reconstructed minus truth on the state's surface. The first two
parameters are named `loc0` and `loc1`, so `writer.histograms()["res_loc0"]` returns
the local-0 residual distribution. Set `geometrySelection` to a list of
`acts.GeometryIdentifier` regions to restrict which states contribute. States
without the requested parameter type are skipped; leaving `parameterType` unset
uses the best available parameters (smoothed, filtered, then predicted).

For residuals against calibrated measurements, set
`reference=acts.examples.TrackParameterReference.Measurement` and explicitly choose
`parameterType`. In that mode none of the truth input collections are required;
only measured coordinates contribute, and efficiency/track-summary plots remain empty.

## Available histograms

The writers return efficiency, fake and duplicate track, track-summary, residual, and pull
histograms according to their configuration. Names and dimensions can vary. To see the available
histograms after `s.run()`:

```python
for name, histogram in sorted(histograms.items()):
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
