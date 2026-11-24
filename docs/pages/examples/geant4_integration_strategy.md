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

## Current integration on this base

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

### Separate session and event ownership

Extract lifecycle and event-result ownership, then introduce configurable action
factories. Retain the serial execution contract and existing public defaults.
Keep the change behavior-preserving and cover both simulation and material
recording. A custom diagnostic action should be installable without editing the
central algorithm, and must not bypass truth or hit bookkeeping.

Acceptance requires unchanged particle and hit output for fixed inputs and
seeds, explicit thread-affinity failures, safe event-result extraction, and
cleanup tests. This is the first implementation step after the baseline.

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

## Decisions for the next design review

Settle the supported initial particle species and generator conventions, the
decay-graph schema, proper-time derivation, and invalid-input policy. Specify
which truth relationships survive filtering and serialization. Choose action
factory ownership and ordering rules that work for serial execution and allow
independent workers later. Select representative benchmark samples and record
the persistent-run baseline before broad performance changes.
