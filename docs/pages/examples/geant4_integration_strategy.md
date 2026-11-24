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

The Examples adapter supports a persistent serial run and an opt-in Geant4
master with native transport workers. Both accept generator-prescribed decay
trees alongside the flat particle collection. See the input contract and
worked configuration below, and the ODD reconstruction measurements at the end
of this page. Material recording continues to use the serial session.

## Changes from the persistent-run baseline

At #6375, `Geant4Simulation` took a flat particle collection and supplied one
Geant4 event per ACTS event. Run initialization, transport and finalization
belonged to the same thread. The run held geometry, materials, regions, cuts
and physics fixed. Random seeds came from the ACTS event and algorithm IDs.

The serial session on this branch preserves those contracts. Both serial and
native sessions now bind the ACTS input to a worker-local event context and
return owned output collections. This keeps the WhiteBoard on the submitting
ACTS thread through output publication. Both modes use the same primary,
tracking, stepping and sensitive-volume mapping code.

The main implementation entry points are listed below.

- [Simulation and action assembly](../../../Examples/Algorithms/Geant4/src/Geant4Simulation.cpp)
- [Persistent run lifecycle](../../../Examples/Algorithms/Geant4/src/RunManager.cpp)
- [Primary translation](../../../Examples/Algorithms/Geant4/src/SimParticleTranslation.cpp)
- [Tracking and truth output](../../../Examples/Algorithms/Geant4/src/ParticleTrackingAction.cpp)
- [Sensitive stepping and hit merging](../../../Examples/Algorithms/Geant4/src/SensitiveSteppingAction.cpp)
- [HepMC input conversion](../../../Examples/Io/HepMC3/src/HepMC3InputConverter.cpp)

## Quasi-stable particle handling

### Configuring a generator-prescribed decay

Use the same HepMC3 conversion for live generators and files. This example
hands charged pion and muon decays to a two-worker Geant4 session. The input
file supplies the decay channels, daughter four-vectors and vertices. Configure
the sequencer with enough ACTS event threads to use the requested workers.

```python
from acts.examples.hepmc3 import HepMC3InputConverter, HepMC3Reader
from acts.examples.simulation import addGeant4

s.addReader(HepMC3Reader(level=level, inputPath=path, outputEvent="gen"))
s.addAlgorithm(HepMC3InputConverter(
    level=level, inputEvent="gen", outputParticles="roots",
    outputVertices="vertices", outputDecayGraph="decays",
    decayPdgCodes=[211, -211, 13, -13],
))
addGeant4(
    s, detector, detector.trackingGeometry(), field, rnd,
    inputParticles="roots", inputDecayGraph="decays", numGeant4Workers=2,
    outputDirRoot=output,
)
```

With `addPythia8`, pass the same `decayPdgCodes` and `outputDecayGraph`
arguments. The generator must be configured to perform the selected decays.
A zero worker count selects serial transport. Omitting both converter options
selects the original flat final-state conversion.

### Input and decay contract

The converter emits selected status-2 parents and the usual unclaimed final
particles in `outputParticles`. The separate `SimParticleDecayGraph` maps a
parent barcode to its proper time, generator endpoint and daughter particles
with their laboratory momenta. All quantities in the graph use ACTS units.
The signed PDG list selects entire decay trees. Each internal species in a
selected subtree must be listed, each vertex must have a single incoming
particle and at least two daughters, and each terminal daughter must have
status 1. Ancestors outside these trees follow the usual HepMC traversal.
Applying a particle selector to the resulting roots removes the whole tree
for each root it rejects. Geant4 requires roots to start inside the
world, or on an inward-facing boundary. An invalid start reports the barcode,
PDG code and position to the submitting event.

For on-shell parents the converter computes proper flight length as
`cΔτ = Δ(ct) m/E`. It checks the generator displacement against
`Δx = Δ(ct) p/E`, checks charge and four-momentum conservation, and rejects
unsupported or inconsistent selected records. Particles need a known PDG
mass and charge, with neutrinos treated as massless. For massive particles,
the converter retains the generator mass when it agrees with the nominal mass
within 5 keV plus 10 parts per million of the nominal mass. This accommodates
differences between particle tables, including the Pythia8 neutral pion.
Off-shell resonances and multi-parent vertices fall outside this contract.
File merging converts each event's space-time coordinates, momenta and
generated masses into the output event's declared units before conversion.

The handover carries decay kinematics. Spin-density matrices and polarization
evolution are outside its scope. The translator attaches the daughters through
`G4PrimaryParticle` and sets the proper decay time. The selected physics
list must supply an active Geant4 decay process for the parent. Transport
determines where decay actually takes
place and boosts daughter momenta using the parent's state at decay. A parent
can bend, lose energy, leave the world, or be removed before decaying.
Daughters enter transport only when Geant4 creates them. Nested prescribed
decays are restored when Geant4 starts tracking the daughter.

### Truth, hits and persistence

A tracked generator daughter carries its original HepMC particle index, an
ACTS barcode reserved for it at conversion, the direct parent barcode and
the Geant4 creation process. The simulated particle starts at the actual
transport vertex. Other Geant4 secondaries use separate barcodes. The merged
HepMC event is the scope of the generator index. Persist the HepMC event
alongside the simulation to keep the original generator vertices.

Hit recording and secondary-killing policies use the Geant4 parent track ID.
A decay inside silicon can produce two distinct hits on that sensor, one for
the incoming parent and one for the outgoing daughter. With
`keepParticlesWithoutHits=False`, transported prescribed parents are retained
so the ancestry of retained daughters can still be followed. The
`addGeant4` ROOT output writes all five components of each direct parent
barcode when a decay graph is supplied. Set
`RootParticleWriter.Config.writeParentBarcodes` explicitly when constructing a writer
by hand. `RootParticleReader` also reads older files, which carry an unknown
direct parent. CSV particle output does not carry these direct-parent fields.

## Adapter responsibilities

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
Start each sample in a fresh process and give it a separate output directory.

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
per serial session or native worker. The built-in kill and hit actions run
first, followed by the configured actions in order. The
[telescope test](../../../Tests/UnitTests/Examples/Algorithms/Geant4/Geant4SimulationActionsTests.cpp)
exercises this interface over two events. It suspends a track before the
first sensor and verifies that its final position and hit count include
the rest of its flight after Geant4 resumes tracking.

The refactor was checked against `e0e07d517d` across the eight benchmark
configurations above. All particle and hit ROOT branches matched bit for bit
over 160 events. Rerunning telescope events 12–14 separately reproduced all
384 hits and the particle records for those events. The material-recording
test reproduced all 60 material tracks both across builds and across event
partitions. The lifecycle, event-context, session-reuse, ordered-action and
sensitive-surface-mapping C++ tests passed on Linux.

### Validate the complete quasi-stable workflow

The [HepMC3 conversion tests](../../../Tests/UnitTests/Examples/Io/HepMC3/HepMC3DecayConversionTests.cpp)
exercise a neutral kaon with two decaying daughters, mixed root and final-state
inputs, alternate units, neutrinos and malformed decay records.
The [end-to-end telescope test](../../../Python/Examples/tests/test_geant4_hepmc_decays.py)
checks nested pion and muon decays across three silicon planes. It compares
serial and native ROOT particle and hit records by event and reads the direct
parent barcodes back through the ACTS reader. Cases include decay inside a
plane, escape before decay, a neutral kaon at rest, and nuclear capture of a
stopped negative pion. A transverse 5 T run checks the daughter vertices,
proper flight times and daughter momentum sums against the bent parent.

A separate 30-event Pythia8 sample used hadronic Z decays at 91.2 GeV,
selected `K_S → π⁺π⁻`, and a 100 mm radial and longitudinal production-vertex
cut in the telescope. Live Pythia and HepMC-file replay both transported 35
selected kaons and realized 29 prescribed decays with the same generator
parent and daughter identities. Serial and two-worker transport of the same
file reproduced all 36 particle and 25 hit branches bit for bit for 1,293
particles and 696 hits. ACTS seeds include the algorithm index, so compare
transport with the same collection and algorithm ordering when changing
worker counts.

In a live Pythia8 top-pair event with the full ODD geometry, Geant4
transported 60 selected neutral pions and produced 120 generator daughters
with direct parent links. The event produced 2,653 hits and 64 CKF tracks. A two-event
`K_S` file replay through the full ODD reconstruction chain yielded 12,124
particles, 705 hits and 26 tracks. Serial and two-worker results agreed bit
for bit across all 160 ROOT fields.

```sh
ctest --test-dir build --output-on-failure -R 'Geant4|HepMC3Decay'
python -m pytest -q Python/Examples/tests/test_geant4_hepmc_decays.py
```

#### Charged-pion decays of the short-lived neutral kaon

In 1,024 Pythia8 hadronic Z events at 91.2 GeV, 1,202 kaons were generated
with the forced decay `K_S → π⁺π⁻`. The ODD input selection retained 1,175.
Geant4 produced 1,136 prescribed pairs, 22 parents reached the transport
volume limit, and 17 produced nuclear secondaries without the pion pair.

The default ODD reconstruction gave unique, charge-correct CKF matches for
both daughters of 102 decays. Truth identified the pairs. The decay-position
estimate uses the closest approach of the two fitted helices in the 2 T
field. Assigning pion masses to their momenta at closest approach gave a
median dipion mass of 497.81 MeV/c², with the central 68% between 495.07 and
500.52 MeV/c². The median 3D distance to the simulated decay vertex was
0.14 mm. All matched pairs enter the plot without a mass or radius cut.

![Dipion mass from the two fitted tracks and reconstructed decay radius](figures/ks-reconstruction.png)

An independent gun sample contained 4,000 kaons at each of 0.5, 2 and
5 GeV/c. Pythia8 sampled the lifetime and generated the decay in each
event. Every kaon decayed inside the zero-field vacuum telescope. The four-worker
and serial outputs had identical ROOT content hashes for particles and hits.
At daughter creation, the mean pion rest-frame momentum was 205.9712 MeV/c.
Two-body kinematics gives 205.9716 MeV/c using the PDG masses. Pythia used
497.610 MeV/c² for the kaon and 139.570 MeV/c² for the pion. The small forward
excess (6,102 of 12,000 positive pions) is already present in the HepMC sample.

![Pion momenta, invariant-mass difference and decay angle at creation](figures/ks-kinematics.png)

The gun's mean proper flight length, measured from the Geant4 vertices and
kaon boost, was 26.83 ± 0.24 mm (statistical). The PDG value is 26.844 mm.
In the ODD sample, five of the 1,136 prescribed decays moved by 0.81–3.94 mm
relative to the generator endpoint. Four parents produced nuclear recoil
tracks and a Geant4 step trace for the fifth recorded elastic scattering in
a cooling pipe. The dipion invariant mass was preserved for all five.
[PDG neutral kaon properties](https://pdg.lbl.gov/2025/tables/rpp2025-tab-mesons-strange.pdf)
and [charged pion properties](https://pdg.lbl.gov/2025/tables/rpp2025-tab-mesons-light.pdf).

![Proper-flight-length distribution and vertex shifts in detector material](figures/ks-flight.png)

To reproduce the three figures, use an ACTS Examples build with ODD, ROOT,
Pythia8, HepMC3 and multithreaded Geant4 with its physics datasets.
The Python analysis uses NumPy, uproot, awkward, pyhepmc and Matplotlib.
The particle gun accesses the installed Pythia8 headers and library through
ROOT's cppyy.

```sh
s=Examples/Scripts/Python
d=build/ks-check
python "$s/geant4_ks_samples.py" generate-gun "$d/gun.hepmc3"
python "$s/geant4_ks_samples.py" transport "$d/gun.hepmc3" "$d/gun" --workers 4
python "$s/geant4_ks_samples.py" generate-z "$d/z.hepmc3" --events 1024
python "$s/full_chain_odd.py" --hepmc3 "$d/z.hepmc3" --decay-pdg-codes 310 \
  --geant4 --geant4-workers 3 -j 4 -n 1024 --output-root -o "$d/odd"
python "$s/geant4_ks_candidates.py" --hepmc "$d/gun.hepmc3" \
  --particles "$d/gun/particles_simulation.root" --output "$d/gun.csv"
python "$s/geant4_ks_candidates.py" --hepmc "$d/z.hepmc3" \
  --particles "$d/odd/particles_simulation.root" --output "$d/odd.csv"
python "$s/geant4_ks_tracks.py" --candidates "$d/odd.csv" --odd "$d/odd" --output "$d/fits.csv"
python "$s/plot_geant4_ks_decays.py" --gun "$d/gun.csv" --odd "$d/odd.csv" \
  --fits "$d/fits.csv" --output "$d/figures"
```

### Improve measured overhead and event throughput

Incorporate #6407 and compare serial and multiprocessing results and resource
use. Preserve global event IDs and seeds, deterministic output ordering, and
worker failure propagation. Test the decay workflow through the same path.

Profile the field cache created on every lookup, the copied input particle
collection, event-store allocation churn, and per-step volume-name checks and
surface lookups. These are optimization candidates, not established dominant
costs. Reuse caches and storage only with explicit context and lifetime rules;
replace mapping or hit-dispatch mechanisms only when profiles justify it.

### Native Geant4 workers and ACTS reconstruction

Set `numGeant4Workers=N` in `addGeant4` or use `--geant4-workers N` in
`full_chain_odd.py` to run Geant4 transport on N native workers. The default
of zero uses the serial session. Native transport requires a multithreaded
Geant4 build. Start a fresh process to create a second native pool, switch
backends or change the physics list. Compatible serial algorithms can reuse
the same run-manager handle sequentially on its owner thread. The ACTS
sequencer thread count controls how many events can run concurrently.

For an input file containing generator-prescribed `K0S` decays, run the
full ODD chain with two native workers and write the parent barcodes to ROOT.
The second command runs the same chain with a live Pythia8 top-pair event,
zero pileup and prescribed neutral pion decays.

```sh
python Examples/Scripts/Python/full_chain_odd.py --hepmc3 events.hepmc3 \
  --decay-pdg-codes 310 --geant4 --geant4-workers 2 -j 4 -n 2 --output-root -o build/odd-ks
python Examples/Scripts/Python/full_chain_odd.py --ttbar --ttbar-pu 0 \
  --decay-pdg-codes 111 --geant4 --geant4-workers 2 -j 4 -n 1 --output-root -o build/odd-pi0
```

The master constructs the detector and physics list and maps sensitive
volumes once. Each native worker owns its actions, field wrapper and active
event state. It obtains the ACTS event number and seed when Geant4 requests
work in `SetUpAnEvent`. The originating ACTS task waits for its own particle
and hit collections. While ACTS digitizes hits or fits tracks for a completed
event, workers can transport other events. At most N requests are admitted
from submission through delivery of the result. The synchronous wait occupies
an ACTS thread, so the number of concurrent ACTS events also limits the
available transport work. See the
[native simulation session](../../../Examples/Algorithms/Geant4/src/NativeSimulationSession.cpp).

Closing the input lets workers drain admitted events and finish the Geant4
run, including when ACTS submitted zero events. The total number of ACTS
events reaching transport can be unknown. Geant4's planned event count for
this streaming run is the worker count used to start the run. Use ACTS event
IDs for event identity. An ACTS action exception reaches the submitting task
after the worker has released its references to the event. The
[session tests](../../../Tests/UnitTests/Examples/Algorithms/Geant4/NativeSimulationSessionTests.cpp)
cover an empty run, startup failure, overlapping events and shutdown. At
shutdown the master joins the workers, detaches detector-region roots while
their logical volumes are live, and releases the geometry on the master thread.

The earlier [handover experiment](../../../Tests/UnitTests/Examples/Algorithms/Geant4/NativeWorkerTests.cpp)
used geantinos in a two-plane telescope. It held one particle inside silicon
until ACTS fitted a straight line for another event, and exercised worker
shutdown with fewer inputs than requested. Its empty terminal slots were
specific to that fixture's handover in `GeneratePrimaries`. The simulation
session ends the worker loop at the event request boundary.

```sh
cmake --build build --target ActsUnitTestGeant4NativeWorkers ActsUnitTestGeant4NativeSimulationSession
ctest --test-dir build --output-on-failure -R '^Geant4Native'
```

#### ODD simulation and reconstruction measurements

The sample used the ODD v6.0.3 tracker, mapped 2 T solenoid field,
`FTFP_BERT` and a 1 m step cap. Each of the 96 events (IDs 5–100) contained
32 muons at 5 GeV with pseudorapidity between 1 and 1.5, using seed 42.
ACTS digitized hits, seeded tracks from space-point triplets, ran the
combinatorial Kalman filter (CKF) and matched tracks to truth. ROOT writers
recorded particles, hits and track summaries including fitted covariances.
A separate debugger snapshot caught the ACTS CKF executing while both native
workers were stepping tracks in Geant4.

| Transport | ACTS threads | Geant4 workers | Events/s | Peak RSS (MiB) | Setup + run (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Serial baseline | 1 | 0 | 9.56 | 1073 | 21.28 |
| Native | 4 | 1 | 10.63 | 1115 | 20.47 |
| Native | 4 | 2 | 20.19 | 1136 | 11.84 |
| Native | 4 | 4 | 30.35 | 1175 | 8.51 |

The table gives medians over three fresh processes per configuration on a
Ryzen 7 5800X with Geant4 11.4.3. Throughput covers event processing and
end-run. Setup plus run also includes detector setup, run initialization and
closing the ROOT files. The serial row used the frozen persistent-run build.
Compared by event ID, every native run reproduced all 155 particle, hit
and track-summary branches bit for bit, including within-event ordering, for
772,019 particles, 48,929 hits and 4,520 tracks. These timings characterize
this muon sample and detector configuration.

From the repository root with the ACTS Python bindings, ODD, ROOT and Geant4
datasets available, run the driver once per configuration in a fresh process.
Set `threads=1` and `workers=0` for serial, or `threads=4` and `workers` to
1, 2 or 4 for native transport. Each run requires a new output directory.

```sh
threads=4
workers=2
python Examples/Scripts/Python/geant4_persistent_run_benchmark.py "build/odd-ckf-${threads}-${workers}" \
  --geometry odd --field mapped --reconstruct --particles 32 --events 96 --skip 5 --seed 42 \
  --acts-threads "${threads}" --g4-workers "${workers}"
```

Use the ROOT content hasher to compare runs on the same build, allowing
events to be written in a different order.

```python
from pathlib import Path
from acts.examples.root import hashRootFile

serial = Path("build/odd-ckf-1-0")
native = Path("build/odd-ckf-4-4")
for name in ("particles_simulation.root", "hits.root", "tracksummary_ckf.root"):
    assert hashRootFile(serial / name) == hashRootFile(native / name), name
```

## Further validation

Extend the generator sample to heavy flavour and high pileup before setting
broader throughput targets. Additional generator conventions and off-shell
parents require a definition of mass, proper time and reference kinematics
for the handover.
