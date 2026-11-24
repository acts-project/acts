@page geant4_integration_strategy Geant4 integration strategy

# Geant4 integration strategy

Evolve the ACTS Examples Geant4 integration around a persistent simulation
session, explicit event ownership, configurable actions, and preservation of
generator decay trees. The first physics milestone is a particle that traverses
detector material and leaves hits before undergoing its generator-prescribed
decay, with its daughters subsequently transported by Geant4.

## Foundation and scope

This development branch starts from
[PR #6375](https://github.com/acts-project/acts/pull/6375), at commit
`de8258c8597e166876e6f3cd97f9d238f8cf68c0`. Its persistent-run lifecycle is the
foundation for subsequent changes. As checked on 24 November 2026, both #6375
and [PR #6407](https://github.com/acts-project/acts/pull/6407) remain open.

#6375 removes repeated run initialization between ACTS events and establishes
explicit begin-run, process-event, and end-run operations. #6407 complements it
with independent spawned processes and output merging. The former addresses
per-event overhead and lifecycle; the latter increases event throughput while
requiring a geometry copy in each worker. Incorporate #6407 after it lands,
without making the decay work depend on multiprocessing.

Keep this strategy and subsequent implementation commits above #6375. If its
head changes, rebase the follow-up work onto the updated head. Once it merges,
rebase only the follow-up commits onto upstream main so the PR's commits are not
duplicated, particularly if it is squash-merged.

The scope is the Examples simulation adapter, its input conversion, and truth
and hit output. Existing geometry conversion and detector construction remain
building blocks. Native Geant4 multithreading, a wholesale simulation rewrite,
and a general external-decayer interface are deferred. The proposed interfaces
below describe responsibilities; their final C++ APIs remain to be designed.

## Integration at the persistent-run baseline

`Geant4Simulation` converts a flat `SimParticleContainer` into G4 primaries,
executes transport, and publishes simulated particles and hits. One ACTS input
event supplies one G4 event containing its particles. The process-wide manager
provides a single serial session. A mutex serializes transport access; it does
not turn the serial run manager into a multithreaded simulation.

On this base, `RunManager` keeps the run open across events. Initialization,
event processing, and finalization must use the same thread. Geometry,
materials, regions, cuts, and physics remain fixed during the run. Per-event
seeds are assigned before transport. These constraints are part of the session
contract and must survive refactoring.

Useful extension points already exist: detector construction, a physics-list
factory registry, an ACTS field adapter, surface mapping, and an internal
stepping-action list. Tracking records parent barcodes and creation processes
for G4 secondaries. However, the simulation constructor installs concrete
primary-generation, tracking, and stepping actions, and these share mutable
event state containing framework data handles and a WhiteBoard pointer.

Output extraction currently occurs after the base execute method releases its
mutex. The supported single-thread contract prevents concurrent event use;
future execution models must make ownership cover output extraction as well as
transport, rather than relying on that mutex alone.

The main implementation entry points are:

- [Simulation and action assembly](../../../Examples/Algorithms/Geant4/src/Geant4Simulation.cpp)
- [Persistent run lifecycle](../../../Examples/Algorithms/Geant4/src/RunManager.cpp)
- [Primary translation](../../../Examples/Algorithms/Geant4/src/SimParticleTranslation.cpp)
- [Tracking and truth output](../../../Examples/Algorithms/Geant4/src/ParticleTrackingAction.cpp)
- [Sensitive stepping and hit merging](../../../Examples/Algorithms/Geant4/src/SensitiveSteppingAction.cpp)
- [HepMC input conversion](../../../Examples/Io/HepMC3/src/HepMC3InputConverter.cpp)

## Quasi-stable particle handling

Geant4 can already decay transported particles through its configured physics.
The missing capability is preserving a generator-prescribed decay while
transporting the parent before that decay.

The current HepMC converter traverses particles with end vertices but retains
only their final-state descendants as simulation particles. It can retain
displaced daughter production vertices, but the parent does not traverse the
intervening detector. The G4 translator then inserts each input particle as an
independent primary, without daughter relationships or prescribed decay times.

Geant4 provides the required transport mechanism: `G4PrimaryParticle` can carry
daughter particles and a prescribed proper decay time through `SetProperTime()`.
Geant4 transports the parent and transforms the prescribed daughters according
to its frame at decay. See the
[event-generator interface](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Fundamentals/eventGenerator.html)
and the
[preassigned decay documentation](https://geant4.web.cern.ch/documentation/pipelines/master/bfad_html/ForApplicationDevelopers/TrackingAndPhysics/physicsProcess.html).

### Input and decay contract

Introduce an optional decay-graph collection alongside the existing particle
collection. Preserve flat inputs for particle guns and existing workflows. The
graph should describe selected transportable parents and their decay trees,
stable generator identities, production information, daughter kinematics, and
prescribed proper decay times. Preserve this information before vertex merging
or final-state selection removes it.

Define selection at the graph boundary: identify which ancestors are
transported, which prompt ancestors are collapsed, and which descendants are
attached to those transported ancestors. An end vertex alone does not establish
a transportable decay; generator records also contain interactions and other
bookkeeping vertices. Validate particle definitions and decay-process support
for every selected species. Exotic particle registration is a separate policy,
not an assumption that every PDG code is supported.

Inject each transport root once. Attach prescribed descendants to their parent;
do not also inject them as independent primaries. Support nested decay chains.
Use proper decay time rather than a fixed laboratory endpoint: energy loss,
scattering, and field bending can change the actual decay position. Define and
validate units, time derivation, rest-frame transformations, and treatment of
inconsistent or incomplete generator records. Keep prescribed decay time
distinct from the particle state's accumulated proper time.

If an interaction destroys a parent before its scheduled decay, its prescribed
daughters must not appear independently later. Preserve original generator
vertices separately from the actual simulated decay vertices.

### Truth and hit contract

Carry generator identity through G4 primary and track information and map actual
G4 track IDs when tracks are created. Preserve identity for prescribed daughters
while allocating separate identities for G4-created particles. Audit the current
translator's preassigned track-ID mapping before extending it to decay trees.

Distinguish generator origin, transport parentage, and creation process. The
current primary/secondary classification uses the presence of a G4 primary
particle pointer, which prescribed daughters can also carry. Use explicit
semantics for secondary killing and hit selection, and retain enough ancestry
when output filters discard particles without hits.

Flush partial sensitive-volume hits when a parent decays or otherwise stops.
The current tracking action warns and clears leftover buffered hits for some
unhandled stop cases. Cover decay, destructive interaction, user termination,
and daughters born inside a sensor; define how each contributes to truth and
hit output.

## Proposed architecture

Keep the ACTS algorithm as a thin adapter with three separable responsibilities.

| Component | Responsibility |
| --- | --- |
| Simulation session | Own geometry, physics, field setup, run lifecycle, and the execution contract |
| Event conversion | Translate flat input or decay graphs, seed the event, and preserve truth identities |
| Actions and output | Collect tracks, hits, material and diagnostics; apply termination policies; return event results |

Construct actions through factories with explicit event context and documented
ordering. Build on the physics-list registry and stepping-action list. Keep
mandatory bookkeeping active when optional actions are added. Factories must
create independent mutable state for each eventual worker; sharing configured
action objects would constrain future concurrency.

Separate framework data access from the state used during transport. A session
should consume one event's input and return its owned result, with no retained
WhiteBoard pointer after completion. Make cleanup and ownership explicit for
event aborts, run aborts, exceptions, and finalization. Material recording must
continue to use the same lifecycle without inheriting unnecessary hit or decay
machinery.

## Delivery sequence and acceptance criteria

### Establish the persistent-run baseline

Complete #6375 and benchmark the resulting integration before attributing
remaining overhead to Geant4 transport. Compare a small detector and a complex
geometry, low and realistic particle multiplicities, and constant and mapped
magnetic fields. Record versions, physics settings, cuts, input, seeds, event
counts, and hardware. Separate startup, run initialization, event conversion,
transport, hit bookkeeping, and output; measure throughput and peak memory.

Acceptance requires reproducible per-event results across event ranges, normal
run/event callbacks, material-output regression coverage, and correct cleanup
after failure. Record numerical performance results before setting optimization
targets.

Use the [serial benchmark](../../../Examples/Scripts/Python/geant4_persistent_run_benchmark.py)
with the ACTS Python bindings, Geant4 datasets and ROOT writers available.
The `odd` geometry also requires the OpenDataDetector and DD4hep bindings.
Start each sample in a fresh process and use separate output directories:

```sh
python Examples/Scripts/Python/geant4_persistent_run_benchmark.py build/sample \
  --geometry telescope --field constant --particles 1 --events 20 --skip 5 --seed 42
```

Repeat with the ODD tracker, mapped solenoid field and higher multiplicity
using `--geometry odd`, `--field mapped` and `--particles 32`. The driver
writes setup and run times and peak process memory to `metrics.json`.
`timing.csv` gives per-algorithm timings for generation, HepMC conversion,
Geant4 simulation and ROOT output. Compare the particle and hit tree contents
for identical event ranges and seeds when changing the adapter.

#### Persistent-run measurements

Measured on 24 November 2026 at `e0e07d517d` with GCC 15.2, a Release build
with assertions, Geant4 11.4.3, HepMC3 3.3.1 and ROOT 6.40.04 on a Ryzen 7
5800X. Each row is one fresh serial process with warm dataset caches, seed 42,
events 5–24, 5 GeV muons with pseudorapidity between 1 and 1.5, FTFP_BERT
default production cuts, a 1 m propagator step cap, and ROOT particle and hit
writers. The telescope has four silicon planes. The ODD v6.0.3 configuration
loads the tracker. Both fields have a nominal central strength of 2 T.

| Geometry | Field | Muons/event | Setup (s) | Initialize (s) | Geant4 algorithm (ms/event) | Event throughput (events/s) | Peak RSS (MiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Telescope | Constant | 1 | 0.37 | 0.34 | 0.079 | 6830 | 433 |
| Telescope | Constant | 32 | 0.37 | 0.32 | 1.270 | 679 | 434 |
| Telescope | Mapped | 1 | 0.43 | 0.33 | 0.115 | 5488 | 434 |
| Telescope | Mapped | 32 | 0.44 | 0.33 | 1.912 | 471 | 435 |
| ODD | Constant | 1 | 1.24 | 9.24 | 2.639 | 361 | 951 |
| ODD | Constant | 32 | 1.37 | 9.36 | 83.308 | 11.5 | 960 |
| ODD | Mapped | 1 | 1.26 | 9.41 | 2.735 | 349 | 951 |
| ODD | Mapped | 32 | 1.26 | 9.30 | 90.429 | 10.5 | 965 |

Setup starts at detector construction. Initialization runs from sequencer
entry through Geant4 run initialization. Throughput covers the subsequent
event window, including generation, conversion, simulation, per-event output
and end-run. Python imports and closing the output files are outside that
window. These short particle-gun runs establish a reference for this host.

For the mapped ODD sample with 32 muons, HepMC conversion took 0.037 ms/event
and the ROOT writers took 4.400 ms/event. A separate instrumented repeat
timed the Geant4 callbacks over the same 20 events. Per event, primary
translation took 0.007 ms, the sensitive stepping action 1.317 ms, the kill
action 1.505 ms and the tracking actions 5.296 ms. The remaining Geant4 event
manager time, including transport and its internal bookkeeping, was 83.007 ms.
The throughput table uses the uninstrumented runs.

### Separate session and event ownership

Extract lifecycle and event-result ownership, then introduce configurable action
factories. Retain the serial execution contract and existing public defaults.
Keep the change behavior-preserving and cover both simulation and material
recording. A custom diagnostic action should be installable without editing the
central algorithm, and must not bypass truth or hit bookkeeping.

Acceptance requires unchanged particle and hit output for fixed inputs and
seeds, explicit thread-affinity failures, safe event-result extraction, and
cleanup tests. This is the first implementation step after the baseline.

The serial implementation reads the particle collection in the ACTS algorithm.
Its [simulation session](../../../Examples/Algorithms/Geant4/src/SimulationSession.hpp)
binds those particles and the geometry context for one synchronous Geant4
event. When transport finishes, the algorithm receives owned result
collections and publishes them to the WhiteBoard. An exception ends the
binding before the event data goes out of scope. Material recording uses the
same session and publishes its material tracks.

For a C++ stepping diagnostic, append a factory to
`Geant4Simulation::Config::steppingActionFactories`. Each factory receives
read-only access to the active event during callbacks and creates one action
per simulation instance. The built-in kill and hit actions run first, followed
by the configured actions in order. The
[telescope test](../../../Tests/UnitTests/Examples/Algorithms/Geant4/Geant4SimulationActionsTests.cpp)
exercises this interface over two events.

The refactor was checked against `e0e07d517d` across the eight benchmark
configurations above. All particle and hit ROOT branches matched bit for bit
over 160 events. Rerunning telescope events 12–14 separately reproduced all
384 hits and the particle records for those events. The material-recording
test reproduced all 60 material tracks both across builds and across event
partitions. The lifecycle, event-context, session-reuse, ordered-action and
sensitive-surface-mapping C++ tests passed on Linux.

### Implement one complete quasi-stable workflow

Start with HepMC input and a controlled parent that crosses sensitive material,
bends and loses energy, decays into prescribed daughters, and leaves connected
truth and hit records. Agree on the decay-graph schema and filtering semantics
before changing the converter. Extend the G4 translator and truth mapping, then
exercise the entire input-to-output path.

Acceptance requires:

- Exactly one transported copy of each intended particle and no duplicate daughters.
- Parent hits before decay and daughter hits after decay, including decay inside a sensor.
- Correct proper decay time and energy-momentum consistency at the simulated decay.
- Preserved identities and ancestry through nested decays and output serialization.
- Neutral parents, decay at rest, and parents leaving the world before decay.
- No prescribed daughters after a destructive interaction removes the parent.
- Defined behavior for missing definitions, invalid decay records, selection, and kill policies.
- Unchanged flat-input behavior and reproducibility across event partitioning.

Reconcile the detailed ATLAS quasi-stable implementation with these contracts
during design review, especially selection, proper-time derivation, and truth
identity. The implementation should use stock Geant4's decay support.

### Improve measured overhead and event throughput

Incorporate #6407 and compare serial and multiprocessing results and resource
use. Preserve global event IDs and seeds, deterministic output ordering, and
worker failure propagation. Test the decay workflow through the same path.

Profile the field cache created on every lookup, the copied input particle
collection, event-store allocation churn, and per-step volume-name checks and
surface lookups. These are optimization candidates, not established dominant
costs. Reuse caches and storage only with explicit context and lifetime rules;
replace mapping or hit-dispatch mechanisms only when profiles justify it.

### Reassess native multithreading

Consider a G4 master/worker backend if multiprocessing's geometry memory or
startup cost becomes limiting. It needs worker-local actions and event state,
deterministic event assignment, and coordination with ACTS scheduling. Compare
throughput, memory, and reproducibility against multiprocessing before adopting
it. Quasi-stable support must not depend on this decision.

#### Native worker handover experiment

The [native worker test](../../../Tests/UnitTests/Examples/Algorithms/Geant4/NativeWorkerTests.cpp)
connects the ACTS TBB sequencer to a stock `G4MTRunManager`. Each ACTS event
submits its ID, seed and input to an exchange and waits for its own future.
The exchange admits at most two queued or running requests. Each Geant4 worker
takes the next ready request in `GeneratePrimaries`, keeps the ACTS identity
in the Geant4 event, and returns owned silicon crossing positions in
`EndOfEventAction`. Reconstruction fits a straight line through those two
positions on the originating event's WhiteBoard.

The two-worker case holds event 17 inside its first silicon step until ACTS
has reconstructed another event. The one- and two-worker runs check the
same per-event random draws and track parameters. A separate run requests
64 Geant4 events and submits four ACTS events. A stepping failure reaches
the originating caller, which skips reconstruction for that event. Closing
the input wakes waiting workers. Each one soft-aborts its own run after at
most one empty terminal event. With Geant4 11.4.3, three ACTS events were
reconstructed and Geant4 counted six slots, including two empty ones.

The following commands run the experiment in a Geant4 MT build with the
example unit tests enabled.

```sh
cmake --build build --target ActsUnitTestGeant4NativeWorkers
ctest --test-dir build --output-on-failure -R '^Geant4NativeWorker'
```

The fixture uses geantinos in a field-free telescope. A production backend
would also need worker-local fields and sensitive actions, validation with
full physics and detector conditions, and a policy for run-size estimates
and failures. The synchronous wait occupies an ACTS event-processing thread.
A full-chain throughput measurement would quantify that cost.

## Decisions for the next design review

Settle the supported initial particle species and generator conventions, the
decay-graph schema, proper-time derivation, and invalid-input policy. Specify
which truth relationships survive filtering and serialization. Choose action
factory ownership and ordering rules that work for serial execution and allow
independent workers later. Select generator and pileup samples to complement
the particle-gun baseline before setting broader performance targets.
