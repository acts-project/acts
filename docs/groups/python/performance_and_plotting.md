@defgroup python_performance_plotting Performance Evaluation and Plotting
@ingroup python_bindings
@brief Extracting efficiency/resolution numbers from a `Sequencer` run, and plotting them.

# ROOT-free performance writers

`acts.examples.PythonPatternRecognitionPerformanceWriter` and
`acts.examples.PythonTrackParameterPerformanceWriter` compute the same efficiency/fake-rate and
residual/pull histograms as their ROOT-backed counterparts (`RootPatternRecognitionPerformanceWriter`,
`RootTrackParameterPerformanceWriter`), but keep the result in memory instead of writing a ROOT
file — the only option available in the PyPI wheel, and often more convenient even with ROOT
available.

```python
import acts.examples.scipy as acts_scipy

cfg = acts.examples.PythonPatternRecognitionPerformanceWriter.Config()
cfg.inputTracks = "fitted_tracks"
cfg.inputParticles = "particles"
cfg.inputTrackParticleMatching = "track_particle_matching"
cfg.inputParticleTrackMatching = "particle_track_matching"
cfg.inputParticleMeasurementsMap = "particle_measurements_map"
writer = acts.examples.PythonPatternRecognitionPerformanceWriter(cfg, acts.logging.INFO)
s.addWriter(writer)

# PythonTrackParameterPerformanceWriter additionally needs a fit backend for the
# residual/pull Gaussian fits -- there is no ROOT-free default, so it must be set explicitly:
cfg_fit = acts.examples.PythonTrackParameterPerformanceWriter.Config()
cfg_fit.inputTracks = "fitted_tracks"
cfg_fit.inputParticles = "particles"
cfg_fit.inputTrackParticleMatching = "track_particle_matching"
cfg_fit.fitFunction = acts_scipy.makeScipyHistogramFitFunction()
fit_writer = acts.examples.PythonTrackParameterPerformanceWriter(cfg_fit, acts.logging.INFO)
s.addWriter(fit_writer)

s.run()

histograms = writer.histograms()  # dict[str, Histogram1 | ProfileHistogram1 | Efficiency1]
efficiency = sum(histograms["trackeff_vs_pT"].accepted.values()) / sum(
    histograms["trackeff_vs_pT"].total.values()
)
```

`Examples/Scripts/Python/pypi_finding_fitting_demo.py` builds and runs both writers end to end
against a custom Python track finder/fitter — the reference to copy from.

`TrackTruthMatcher(doubleMatching=True)` is the standard way to produce the
`inputTrackParticleMatching`/`inputParticleTrackMatching` collections these writers need.

> [!warning]
> Both writers evaluate over the *unfiltered* `inputParticles` collection, not just matched
> tracks — pass an already-selected particle collection (e.g. `particles_selected`), or the
> writer step can dominate your total runtime on realistic event sizes. See @ref python_caveats.

# Plotting

`Histogram1`/`Histogram2`/`Histogram3`, `ProfileHistogram1`, and `Efficiency1` — the objects
`.histograms()` returns — support `.plot(ax=None, **kwargs)` directly (matplotlib + mplhep):

@snippet{trimleft} examples/test_performance_and_plotting.py Plotting an ACTS histogram

They also convert to a [boost-histogram](https://boost-histogram.readthedocs.io/) object, so any
tool that consumes one works unmodified (rebinning, further `mplhep` styling, `hist`'s plotting
shortcuts, ...):

@snippet{trimleft} examples/test_performance_and_plotting.py Converting to boost-histogram

For geometry and track visualization (not histogram plotting), see
`acts.examples.visualization.PyVisualization2D` and `TrackVisualizerAlg`, and
`Examples/Scripts/generic_plotter.py` for the YAML-configured plotting tool used by physmon.

> [!note]
> `matplotlib`, `mplhep`, and `boost_histogram` are optional dependencies for this — install
> them alongside `pyacts` if you want to plot.
