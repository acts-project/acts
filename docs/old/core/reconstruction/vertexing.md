(vertexing_core)=
# Vertexing

:::{todo}
Write vertexing documentation
:::

## Displaced two-track prototype

`Acts::estimateTwoTrackVertexSeed` estimates a spatial seed from two tracks in a
uniform axial magnetic field. It uses transverse circle intersections and the
longitudinal separation of the tracks to resolve the two-intersection ambiguity.
For nonintersecting circles it uses their closest transverse points; in zero
field it uses closest points on straight lines. Parallel lines and concentric
circles without a unique seed return a seeding error. The result contains a
position and a three-dimensional point separation, not a covariance or an exact
helix-to-helix distance of closest approach.

The circle-based approach follows the geometric strategy of Athena's
`Tracking/TrkVertexFitter/TrkVKalVrtCore/src/VKvFast.cxx` (`vkvFastV`), implemented
with ACTS charge, field and parameter conventions. It does not include the
reference-primary-vertex pointing heuristic. Only the nearest half-turn of each
trajectory is considered; looping tracks and ambiguous crossings can need a
more complete seeder.

`Acts::TwoTrackVertexFinder` composes this estimator with `Acts::IVertexFitter`.
It implements `Acts::IVertexFinder` and fits each pair independently:

```cpp
Acts::TwoTrackVertexFinder::Config cfg;
cfg.vertexFitter = fitter;  // shared_ptr<const IVertexFitter>
cfg.field = fieldProvider;  // shared_ptr<const MagneticFieldProvider>
cfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
Acts::TwoTrackVertexFinder finder(cfg);
auto state = finder.makeState(mctx);
Acts::VertexingOptions options(gctx, mctx);  // Unconstrained fit.
auto result = finder.find(inputTracks, options, state);
```

The default selection requires opposite charges, seed separation at most 1 mm,
and total fit chi-squared at most 25. These are prototype cuts, not tuned physics
working points. The finder samples the field at the midpoint of the input track
positions and approximates it as uniform and axial for seeding. It rejects
non-axial fields; the fitter remains responsible for its own propagation. Its
propagator and the finder's field provider must describe the same field.

Keep the field context alive and unchanged while using the state, and use the
same context in `VertexingOptions`. Construct a fresh state for each event.
`setTracksToRemove` excludes tracks on subsequent calls with that state. Inputs
are deduplicated by track identity; output vertices may still share tracks.
Seed degeneracies reject the pair. Other seed errors, field errors and fitter
errors fail the call rather than silently hiding failures.

There is no primary-vertex selection, mass constraint, mass tagging, overlap
resolution or vertex merging. In particular, do not pass a primary vertex as a
fit constraint: constrained `VertexingOptions` are rejected. Runtime scales
quadratically with input multiplicity. This is a small reconstruction building
block, not a production secondary-vertex selection.

### Synthetic displaced pion validation

The `TwoTrackVertexFinder` unit-test executable includes a deterministic synthetic
particle gun. It creates back-to-back pion daughters in a 497.611 MeV parent's
rest frame, boosts them with parent momenta of 1, 2 and 5 GeV, and places their
common production vertex at transverse radii of 1, 10, 50 and 100 mm. The parent
moves transversely; vertex z is sampled within ±20 mm. The decay angular sample
is restricted to rest-frame |cos(theta)| < 0.8.

The daughters are propagated with an Eigen stepper in a uniform 2 T axial field
to a common origin perigee. Before smearing, the geometric seed must recover the
known production point. The perigee measurements then receive independent
Gaussian errors: 20 µm in d0, 50 µm in z0, 0.2 mrad in each angle, 0.1% in q/p,
and 0.1 ns in time. These assumed diagonal covariances are also given to Billoir.
The finder receives neither the true vertex nor a mass hypothesis constraint.

Run the scan and print its per-bin acceptance and aggregate residual, pull and
mass diagnostics with:

```sh
cmake --build build --target ActsUnitTestTwoTrackVertexFinder
build/bin/ActsUnitTestTwoTrackVertexFinder --log_level=message
```

The tests check analytic geometry cases, interface dispatch and error forwarding,
track exclusion, simple unrelated-track rejection, and 240 smeared signal pairs.
They require at least 16 accepted pairs in every 20-event bin, a spatial residual
RMS below 2 mm, a mean mass within 3 MeV of the injected mass, and spatial pull
RMS values between 0.5 and 1.5. These are broad regression checks for this fixed
toy sample, not a calibration of reconstruction uncertainties.

This sample contains no parent transport, decay-radius distribution, detector
hits, material interactions, pattern recognition, or realistic combinatorial
background. Its acceptance measures seed-and-fit behavior for supplied tracks,
not detector efficiency. Realistic K-short simulation and detector-level
validation are separate follow-up work.
