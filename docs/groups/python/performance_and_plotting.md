@defgroup python_performance_plotting Performance Evaluation and Plotting
@ingroup python_bindings
@brief Extracting efficiency/resolution numbers from a `Sequencer` run, and plotting them.

# ROOT-free performance writers

`acts.examples.PythonPatternRecognitionPerformanceWriter` produces efficiency, fake-rate,
duplication, and track-summary histograms. `PythonTrackParameterPerformanceWriter` produces
residual, pull, efficiency, and track-summary histograms. Both expose the results in memory with
`.histograms()`; neither needs ROOT.

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

The track-parameter writer also needs `cfg.fitFunction` for its Gaussian fits; use
`acts.examples.scipy.makeScipyHistogramFitFunction()` without ROOT. See
`Examples/Scripts/Python/pypi_finding_fitting_demo.py` for both writers in a complete chain.

`TrackTruthMatcher(doubleMatching=True)` is the standard way to produce the
`inputTrackParticleMatching`/`inputParticleTrackMatching` collections these writers need.

## Available histograms

The keys depend on the writer and its configuration. This table groups the main families;
`pT`, `eta`, and `phi` variants are often available alongside the examples shown.

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
