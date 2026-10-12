# Calorimeter prototype

Optional calorimeter reconstruction for the `acts-nodd` development line, hosted
in ActsExamples. Detector descriptions and readout encodings belong to the
external [nODD repository](https://github.com/asalzburger/nodd). This package owns
reconstruction helpers, adapters, configuration and validation.

## Executable baseline

Enable `ACTS_BUILD_CALORIMETER=ON` in a normal ACTS build. The option defaults to
OFF and enables the existing basic Examples dependencies. It adds no new
third-party package, source download, detector requirement or Gaudi dependency.
The library targets are `Acts::Calorimeter` (framework-independent reconstruction) and
`Acts::ExamplesCalorimeter` (ActsExamples algorithms).
The `github-ci` preset enables this package and its CTest entries; the existing
coverage preset explicitly disables it to keep its Examples-free configuration.

```sh
cmake -S . -B build/calo -DACTS_BUILD_CALORIMETER=ON \
  -DACTS_BUILD_UNITTESTS=ON
cmake --build build/calo --target ActsExampleCalorimeter \
  ActsUnitTestCalorimeterResponse ActsUnitTestCalorimeterClustering
ctest --test-dir build/calo -R '^Calorimeter(Response|Clustering|Example)$' --output-on-failure
build/calo/bin/ActsExampleCalorimeter
```

Use the same compiler and dependency setup as the rest of ACTS. The example runs
three synthetic events through the real ActsExamples sequencer:

```text
SyntheticDeposits -> CalorimeterDigitizationAlgorithm -> CheckCalibratedHits
                  -> CalorimeterClusteringAlgorithm -> CheckClusters
```

Each event produces calibrated cells 42, 43 and 46 with energies 0.06, 0.04 and
0.06 GeV. The synthetic topology connects cells 42 and 43; cell 46 is isolated.
The two resulting clusters contain 0.10 and 0.06 GeV, conserving the accepted cell
energy. The executable checks both stages and returns failure on a mismatch.
It writes the framework's
timing CSV to the working directory. It neither loads nODD nor simulates showers;
its purpose is to establish a runnable reconstruction and event-store boundary.

## Data and response contract

All values use ACTS native units. Use `Acts::UnitLiterals`/`UnitConstants` when
converting values from another package, including time: do not assume numeric
nanoseconds already have the native ACTS time scale.

- `SimCalorimeterHit` contains a full 64-bit readout cell ID, the global cell
  centre, deposited energy and deposit time. It is a transient reconstruction
  input; it does not replace EDM4hep's persistent simulation data model.
- `CalorimeterHit` contains calibrated energy, an energy-weighted time and
  indices into the input deposit collection. Its position is a cell centre,
  never a shower-step position. Accepted deposits sharing an ID must use exactly
  the same centre, as obtained from one geometry lookup.
- `CalorimeterResponse` validates finite positions/times and finite nonnegative
  deposited energies. It accepts an inclusive deposit-time window, ignores
  zero-energy deposits, merges accepted deposits by cell, applies one positive
  global energy scale, then accepts positive cell energies at or above the
  calibrated threshold. The threshold is applied after merging. Output is
  sorted by cell ID; malformed data or numerical overflow causes failure.
- The time is weighted by accepted deposited energies; neither time of flight
  nor an event-time correction is applied. The scale is an explicit prototype
  parameter, not an established nODD calibration.
- Each event owns its output. The response has immutable configuration and no
  shared event state. Input indices remain meaningful only with the named input
  collection for that event; podio adapters create persistent relations
  explicitly.

## Clustering contract

`CalorimeterClusterer` implements a seeded connected-component baseline. Its
geometry-supplied `CellNeighbour` edges form an undirected graph copied and
normalized once at construction; duplicate/reversed edges are harmless, while
self-edges are rejected. Full readout IDs are preserved. Numeric proximity of IDs
does not imply adjacency. A cell missing from the graph is isolated; a geometry
cell without an accepted event hit cannot bridge components.

The algorithm keeps positive cells at or above the calibrated neighbour-energy
threshold and retains each connected component containing a cell at or above the
seed-energy threshold. The seed threshold must be at least the neighbour
threshold, and both must be finite and nonnegative. A component with multiple
seeds produces one cluster; a component whose total energy exceeds the seed
threshold but has no qualifying individual cell is dropped. Zero-energy and
below-threshold cells do not join or connect clusters.

Each retained cell appears exactly once. Cluster energy is the sum of calibrated
cell energies, with no second response or calibration. Position and time are
calibrated-energy-weighted cell means. The highest-energy cell identifies the
seed, with smaller IDs breaking ties. Clusters are ordered by their smallest
constituent ID; constituent indices and accumulation order follow cell-ID order.
These conventions make results reproducible under hit and edge permutations.

`CalorimeterCluster::hitIndices` refers to the input calibrated-hit collection
for that event. Following a hit's `sourceIndices` reaches simulated deposits.
The ActsExamples adapter writes a new cluster collection and preserves its
input. The native helper has no shared event state. Duplicate input cell IDs,
negative/nonfinite energies, nonfinite positions/times and output overflow are
errors. The IO adapters below translate these indices to persistent relations.

This baseline does not split nearby showers within one connected component,
perform noise-significance clustering, infer geometry, correct time of flight,
or correct shower energies. Layer/readout boundaries must be respected by the
geometry that supplies the edges. The synthetic graph is not a nODD geometry
model. Noise, sampling fluctuations, detector-specific response, saturation,
detailed electronics and shower-energy corrections are later milestones. No jet
or particle-flow performance is established by this synthetic example.

## Optional calorimeter jets

Enable both `ACTS_BUILD_CALORIMETER=ON` and the existing
`ACTS_BUILD_EXAMPLES_FASTJET=ON` option to build `Acts::CalorimeterJets` and
`Acts::ExamplesCalorimeterJets`. This reuses ACTS's existing FastJet discovery;
FastJet must already be installed (ACTS currently requires at least 3.4.1).
The approved local prototype uses FastJet 3.5.1, licensed under
[GNU GPL v2 or later](https://www.fastjet.fr/about.html). No source is copied or
downloaded. With the FastJet option OFF, the response/clustering libraries,
tests and executable continue to build without a FastJet dependency.

```sh
cmake -S . -B build/calo -DACTS_BUILD_CALORIMETER=ON \
  -DACTS_BUILD_EXAMPLES_FASTJET=ON -DACTS_BUILD_UNITTESTS=ON
cmake --build build/calo --target ActsExampleCalorimeterJets \
  ActsUnitTestCalorimeterJets
ctest --test-dir build/calo -R '^Calorimeter(Jets|JetExample)$' --output-on-failure
build/calo/bin/ActsExampleCalorimeterJets
```

`CalorimeterJetReconstruction` assigns each positive-energy cluster the massless
four-momentum `(E * unit(position - origin), E)`, in `(px, py, pz, E)` order and
ACTS native units. The default origin is the global zero point; a configured
origin is fixed for the run. Cluster times are not used. No additional energy
calibration is applied. FastJet inclusive anti-kt uses radius 0.4 by default
in rapidity-phi space and E-scheme four-vector addition. Recombined jets can have
nonzero mass. The configured minimum pT is inclusive and defaults to zero.
Outputs are ordered by decreasing pT, with sorted constituent indices breaking
exact ties for a fixed input collection. FastJet determines clustering in
degenerate configurations; no general permutation-invariance claim is made.

`CalorimeterJet::clusterIndices` refers to the input cluster collection, including
its original indices when zero-energy clusters are skipped. Four-momenta and
constituent indices are copied while the clustering sequence is alive, so event
output holds no FastJet objects or pointers. All clustering state is local to
one call. The ActsExamples algorithm retains its input and writes a new jet
collection. Linking the installed `Acts::CalorimeterJets` component finds
FastJet again through the ACTS package configuration.

Configuration requires a finite positive radius within FastJet's supported range,
a finite nonnegative pT cut, and a finite origin. Input energies must be finite
and nonnegative, and positions finite, including for zero-energy clusters.
Zero-energy clusters are skipped; a positive-energy cluster at the origin has
an undefined direction and is rejected. Undefined/nonfinite displacements,
indices outside FastJet's signed-int range and excessive total energies are
errors. Total energy is bounded by half the square root of the largest double
to keep FastJet's squared-momentum calculations finite. These are numeric
validity checks, not detector acceptance cuts.

`ActsExampleCalorimeterJets` extends the same three-event response and clustering
fixture with `CalorimeterJetAlgorithm -> CheckJets`. The opposing clusters
produce jets with pT 0.10 and 0.06 GeV. Unit tests additionally check nearby
clusters merging, four-momentum conservation and jet mass, the phi wrap,
radius changes, forward-cluster pT ordering, a displaced origin, threshold
boundaries, malformed inputs and the complete deposit-to-jet provenance chain.
This is a calorimeter-only baseline. Track matching, particle flow, pileup
subtraction, detector calibration and realistic jet performance remain later
milestones.

## Optional EDM4hep hit IO

Enable `ACTS_BUILD_CALORIMETER_EDM4HEP=ON` to build
`Acts::ExamplesCalorimeterEDM4hep`. This enables the calorimeter and podio examples
options, which in turn enable the existing
EDM4hep plugin, its data-model dictionary and ROOT IO dependencies. The broader
`ACTS_BUILD_EXAMPLES_EDM4HEP` option is unnecessary for this calorimeter adapter;
no DD4hep or Gaudi adapter is required. The default calorimeter build continues
to work with the IO option OFF. This adapter requires EDM4hep 1.0 or newer and
podio 1.7 or newer for typed calorimeter links; older versions remain supported
by the existing ACTS EDM4hep plugin when the calorimeter IO option is OFF.
Packages must already be installed; no automatic download
or vendoring is added. The user-approved local prototype was tested with
EDM4hep 1.1.1, podio 1.8.0 and
ROOT 6.40.04. EDM4hep and podio use Apache-2.0 licenses
([EDM4hep](https://github.com/key4hep/EDM4hep/blob/main/LICENSE),
[podio](https://github.com/AIDASoft/podio/blob/master/LICENSE)).

```sh
cmake -S . -B build/calo -DACTS_BUILD_CALORIMETER_EDM4HEP=ON \
  -DACTS_BUILD_UNITTESTS=ON
cmake --build build/calo --target ActsExampleCalorimeterEDM4hep \
  ActsUnitTestCalorimeterEDM4hep ActsUnitTestCalorimeterEDM4hepReconstruction
ctest --test-dir build/calo -R '^CalorimeterEDM4hep(Reconstruction|Example)?$' --output-on-failure
build/calo/bin/ActsExampleCalorimeterEDM4hep /tmp/calorimeter-io-example
```

The input converter reads one named `edm4hep::SimCalorimeterHitCollection` from
the podio frame. Each `CaloHitContribution` becomes a native simulated deposit,
in hit order followed by contribution order. GeV, mm and ns are converted
explicitly into ACTS units. Full unsigned 64-bit cell IDs are preserved.
Contribution step positions and simulated-hit positions are not assumed to be
cell centres. A required `cellCentre(cellId)` callback supplies global centres
in ACTS units, with one lookup per unique cell per event. Parallel processing
requires the callback to support concurrent calls. A geometry adapter belongs
to the following nODD milestone.

Input hit energy must match the sum of contribution energies within 16 float
epsilons relative to the larger sum; contributions are never rescaled. A positive
hit without contributions is an error because it provides no deposit times.
Empty collections and zero-energy hits without contributions are accepted.
Negative/nonfinite energies, nonfinite contribution times, invalid/repeated
contribution relations, missing/wrong collections and unresolved/nonfinite
cell centres cause failure. The converter writes a source mapping alongside
the native deposits, retaining the actual EDM4hep hit and contribution handles.

The output converter writes `edm4hep::CalorimeterHitCollection` and
`edm4hep::CaloHitSimCaloHitLinkCollection`. Energies, times and centres are
converted back to GeV, ns and mm, with normal EDM4hep float rounding. Overflow
and positive-energy underflow are errors. No second calibration is applied.
Each output hit requires valid, unique deposit provenance matching its cell.
The source mapping is checked against original contribution values and relations.
One link per contributing simulated hit carries its fraction of the accepted
deposited energy for that calibrated cell; fractions sum to one. Rejected
contributions do not enter weights. These are simulated-hit-level associations,
not separate persistent links to each accepted contribution.

Use the existing `PodioReader` and `PodioWriter`. Configure the writer with the
original `inputFrame` and both names from the output converter's `collections()`.
Retaining the frame preserves original simulated hits, contributions and MC
particles so all persistent relations resolve after writing and rereading.
Output collection names must also be absent from the input frame. Native
deposits, calibrated hits and source mappings remain available to subsequent
reconstruction algorithms before the writer consumes the frame.

The example writes `calorimeter-input.root` and `calorimeter-output.root` in its
supplied directory and checks a three-event ROOT round trip through the real
ActsExamples reader, converters, response and writer. Its generated fixture
contains a high-bit cell ID, two simulated hits in the same cell, contribution
times, a rejected late contribution and MC relations. It verifies calibrated
energies, time, geometry centres and truth weights 0.75/0.25 after rereading.
The extended fixture adds two cells: one joins the high-ID cell in a 0.12 GeV
cluster, and one makes a separate nearby 0.06 GeV cluster. The opposing cluster
has 0.04 GeV. It also checks persistent cluster relations, times, seed IDs and
reconstruction metadata. This is a synthetic EDM4hep fixture, not a detector
simulation sample.

## Persistent clusters, jets and configuration

`EDM4hepCalorimeterClusterOutputConverter` reads native clusters, native hits and
the EDM4hep hits from the hit converter. It preserves native collection order
and writes `edm4hep::ClusterCollection` with `Cluster::hits` relations. Input
hit mapping is checked by cell ID and converted energy, time and position;
indices alone are insufficient to identify a persistent object. Cluster energy,
centre, time and seed must agree with its unique constituent hits. A hit cannot
belong to two clusters. Invalid indices, duplicate cells, inconsistent values
and float overflow cause failure before any output collection is stored.

EDM4hep Cluster has no time or seed-ID field. The converter therefore writes a
`podio::UserDataCollection<double>` for times in ns and a
`podio::UserDataCollection<uint64_t>` for full seed IDs. Both have one entry per
cluster in the same order, including empty events. Their names are returned by
`collections()` and recorded in metadata. Never filter or reorder a cluster
collection independently of its sidecars. Intrinsic shower direction, covariance
and shape fields retain their EDM defaults and are marked unmeasured in metadata.

When the IO and FastJet options are both ON, the additional
`Acts::ExamplesCalorimeterJetsEDM4hep` library provides
`EDM4hepCalorimeterJetOutputConverter`. It writes a named
`edm4hep::ReconstructedParticleCollection` for jets with energy, momentum and mass
in GeV and `ReconstructedParticle::clusters` relations. The collection name
identifies the objects as jets; PDG, charge and covariance are not measured.
Momentum and energy are rounded to EDM4hep floats separately. Mass is computed
from the native four-vector, allowing roundoff-sized negative mass squared to
clamp to zero. No jet constituent particles or PFOs are invented.

Use the same jet reconstruction configuration for the native jet algorithm and
its converter. The converter verifies native/persistent cluster mappings and
E-scheme sums using that fixed origin, unique cluster membership and the pT cut.
A small double-precision tolerance scales with constituent count to allow
FastJet's different summation order. Native constituent order and output jet
order are retained. The converter adds no second calibration or clustering.

```sh
cmake -S . -B build/calo -DACTS_BUILD_CALORIMETER_EDM4HEP=ON \
  -DACTS_BUILD_EXAMPLES_FASTJET=ON -DACTS_BUILD_UNITTESTS=ON
cmake --build build/calo --target ActsExampleCalorimeterJetsEDM4hep \
  ActsUnitTestCalorimeterEDM4hepJets
ctest --test-dir build/calo -R '^CalorimeterEDM4hep(Jets|JetExample)$' --output-on-failure
build/calo/bin/ActsExampleCalorimeterJetsEDM4hep /tmp/calorimeter-jets-io-example
```

The jet ROOT example merges the two nearby clusters into a 0.18 GeV jet with
nonzero mass and keeps the opposing 0.04 GeV jet. After rereading it checks the
complete jet-to-cluster-to-hit-to-simulated-hit-to-contribution-to-MC chain,
four-momentum and configuration. With FastJet discovery disabled, the hit/cluster
IO library and ROOT example still build and run. No further package is required.

`EDM4hepCalorimeterMetadata` runs after converters and transfers the original
frame from `inputFrame` to a distinct `outputFrame` key. Point PodioWriter's
`inputFrame` to that output key, and include all converter `collections()` names
in the writer. The frame retains simulation collections and existing parameters.
It receives schema version 1 parameters under `acts.calo.`: units, collection
names, response scale/threshold/time window, clustering thresholds and normalized
neighbour edges, geometry identifier and field conventions. Neighbour IDs are
stored as decimal `first:second` strings, preserving all 64 bits. Jet converter
`metadata()` supplies algorithm, radius, pT cut, origin and constituent conventions
through the metadata configuration's `additional` field.

Construct metadata from the same configurations passed to the algorithms and
converters, as the example does; it checks their collection-name connections but
does not inspect the previously executed algorithms. A versioned geometry source
identifier is required because a cell-centre callback cannot be serialized.
Unbounded response time limits are stored as double infinities. Existing keys
under `acts.calo.` and duplicate additional keys across parameter types are
rejected, avoiding stale configuration. Parameters are repeated per event in
this prototype, including topology; a run-level representation for large detector
geometries is a later optimization. No full detector calibration or physics
performance is established by this IO fixture.

## Package layout

```text
Calorimeter/
  include/ActsCalorimeter/   # Framework-independent data and reconstruction
  src/Digitization/         # Deterministic response baseline
  src/Clustering/           # Seeded connected-component baseline
  Examples/                 # ActsExamples algorithms and runnable example
  Jets/                     # Optional FastJet helper and owned jet event data
    Examples/               # ActsExamples jet adapter and extended example
  Io/EDM4hep/               # Optional hit/cluster IO, sidecars and frame metadata
    Jets/                   # Optional jet IO and extended ROOT example
Tests/UnitTests/Calorimeter/ # Response, clusters, jets and integration checks
```

Add geometry adapters, tracking and Pandora subdirectories
when the corresponding implementation lands. Keep dependencies optional at
their point of use; do not introduce empty backends or automatic downloads.

## Reviewable PR sequence

Each stage gets a dedicated branch pushed to `asalzburger/acts`. Stack dependent
PRs explicitly and keep validation evidence in their descriptions. PR 1 starts
from `codex/nodd-tracker-gen3`, whose existing PR provides the nODD integration.

| PR | Scope | Acceptance criterion | External decision |
| --- | --- | --- | --- |
| 1 (merged: #4) | Optional package, data contract, deterministic response, ActsExamples synthetic run | Enabled/disabled build checks, response tests and executable pass | No new packages |
| 2a (merged: #5) | Calorimeter cluster baseline | Seed/neighbour boundaries, energy accounting, provenance, input-order invariance and sequencer checks | No new packages |
| 2b (merged: #6) | Reconstructed calorimeter jets | Geometry/mass convention, FastJet constituent mapping and synthetic jet checks | User approved reuse of installed FastJet 3.5.1 through the existing ACTS integration |
| 3a (merged: #7) | EDM4hep/podio simulated-hit input and calibrated-hit output | ROOT fixture round trip; preserve 64-bit IDs, units and hit/contribution/particle relations | User approved reuse of installed EDM4hep 1.1.1, podio 1.8.0 and ROOT 6.40.04; no vendoring |
| 3b (this change) | Persistent cluster/jet output and reconstruction metadata | Resolve cluster-to-hit and jet-to-cluster relations after rereading; record conventions/configuration | Reuse the approved IO packages |
| 4 | nODD calorimeter geometry adapter and full-simulation input | Cells, centres, layers and neighbours agree with DD4hep; reproducible single-particle samples | Confirm geometry/input source and any needed external acquisition |
| 5 | ACTS track extrapolation to calorimeter entrance surfaces | Barrel/endcap states, covariance, failures and track-hit provenance validated | Reuse ACTS propagation; no Gaudi tracking wrapper |
| 6 | Direct PandoraSDK/LCContent ActsExamples adapter | Geometry/plugins/settings initialization, event conversion, PFO ownership and reliable reset; serial processing first | User approval required before fetching, linking or pinning Pandora dependencies |
| 7 | PFO jets, calibration and physics validation | Compare truth, calorimeter and PFO jets; single-particle response and pileup study with recorded inputs/configuration | Review any further response/clustering packages before adding them |

The first jet baseline can use synthetic calibrated cells; realistic performance
depends on the later geometry and simulation stages. Pandora's lepton-collider
content requires dedicated calibration and validation for proton collisions.

## Reuse policy

Use DD4hep/DDG4, EDM4hep/podio, FastJet and Pandora as library or IO interfaces
where appropriate. Key4hep calorimeter and Pandora packages are reference
implementations for adapters and response algorithms; their Gaudi algorithms
cannot be linked unchanged into this package. Any extracted implementation needs
a review of its license, attribution, framework coupling and detector assumptions.

The user requested a check-in before introducing externals. Researching an
interface or documenting a candidate does not authorize downloading, installing,
vendoring, pinning or linking it. Present the proposed package/version, purpose,
license and build impact for that decision. Existing ACTS dependencies used by
the basic Examples framework are sufficient for PR 1.
