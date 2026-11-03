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

## Regional matching quality

Both pattern-recognition writers can resolve purity and completeness by named
geometry regions. Set `inputMeasurements` to the **complete** measurement collection
indexed by the truth associations and track source links, and configure
`trackQualityGeometrySelections`:

```python
cfg.inputMeasurements = "measurements"
cfg.trackQualityGeometrySelections = {
    "pixel": [acts.GeometryIdentifier(volume=10)],
    "strip": [acts.GeometryIdentifier(volume=20), acts.GeometryIdentifier(volume=30)],
    "pixel_layer": [acts.GeometryIdentifier(volume=10, layer=2)],
}
```

The identifiers above are illustrative; use the identifiers of your geometry.
Selections use the same hierarchical matching as the track-parameter collector's
`geometrySelection`: a volume includes its layers and sensors; a layer includes
its sensors. An empty selection includes all measurements. Regions may overlap;
each is evaluated independently. Region names must be non-empty and prefix the
output names, for example `pixel_completeness_vs_pT` and `pixel_purity_vs_eta`.
The existing `subDetectorTrackSummaryVolumes` remains a separate configuration
for measurement, hole and outlier counts.

These profiles keep the **global** track-to-particle match, including duplicates;
there is no regional rematching or additional quality cut. Fakes and unassociated
tracks do not contribute. For each region:

- Completeness is the number of shared, non-outlier measurements divided by the
  matched particle's measurements in that region.
- Purity uses the same numerator divided by the track's measurements plus
  outliers in that region, consistent with the global quality profiles.

The geometry comes from measurement indices (`IndexSourceLink`), so unfitted
tracks do not need reference surfaces. Holes do not contribute to either ratio.
Profiles use the matched particle's **truth** pT, eta and phi, and the binning of
`trackQualityPlotToolConfig`. If a denominator is zero, that metric is undefined
and is excluded from its profile. A defined ratio with no shared hits is filled
with zero. Completeness and purity can therefore have different entry counts.

Python `.stats()["regionalQuality"][name]` exposes `nCompleteness`, `nPurity`,
`nUndefinedCompleteness`, `nUndefinedPurity`, `meanCompleteness` and `meanPurity`.
An undefined mean is NaN. ROOT writes the same quantities with a `name_` prefix.
An empty region map disables this diagnostic and requires no new input collection;
the global histograms and matching statistics are unchanged.

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
