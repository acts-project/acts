@page examples_geometry_contexts Simulation and reconstruction geometry contexts

# Simulation and reconstruction geometry contexts

Alignment studies need simulation to place the detector modules differently from
what reconstruction assumes. The Examples framework supports this directly: an
`ActsExamples::AlgorithmContext` carries two geometry contexts, and one sequence
can run both geometries at once.

There is no need to write measurements out of a simulation job and read them back
into a separate reconstruction job. Measurement IO stores local parameters and a
geometry identifier, so it carries no geometry context anyway - the two-context
model expresses the same thing without the round trip.

## The two contexts

| Member | Meaning |
| --- | --- |
| `recoGeoContext` | The geometry reconstruction assumes, i.e. the current alignment hypothesis. |
| `simGeoContext` | The geometry the detector is actually built with, i.e. what simulation transports particles through. |

Both default to an empty context. A job that adds no context decorator behaves
exactly as if there were a single context.

The rule for picking one is:

- A quantity derived from **simulation truth** - sim hits, truth particle
  positions on sensitive surfaces - uses `simGeoContext`. Simulation and
  digitization are entirely on this side.
- A quantity derived from **reconstruction output** - track states, fitted
  parameters, space points, measurements as they are placed for pattern
  recognition - uses `recoGeoContext`. Seeding, track finding, fitting,
  extrapolation and vertexing are entirely on this side.
- A chain that is neither, such as material mapping, standalone propagation or a
  geometry dump, uses `recoGeoContext`.

Performance writers legitimately use both. `ActsExamples::RootTrackStatesWriter`
is the clearest case: its truth branch reads sim hits in `simGeoContext` while
its measurement and track state branches sit in `recoGeoContext`, and the
difference between them is the misalignment under study.

Writers that intersect a truth particle with a perigee or beamline surface are
insensitive to the choice, because those surfaces carry no detector element and
therefore no alignment payload.

## Injecting a misalignment

`ActsExamples::AlignmentDecorator` decides which context(s) it writes through
`Config::target`:

| Target | Effect |
| --- | --- |
| `eSim` | Only simulation sees the alignment. Reconstruction stays nominal - this is the alignment study setup. |
| `eReco` | Only reconstruction sees it. |
| `eBoth` | Both, i.e. a detector that is misaligned but perfectly known. This is the default. |

Chain two decorators to give simulation and reconstruction two *different*
non-nominal alignments, which is what an alignment iteration looks like.

## Caveats

- **Geant4 cannot be misaligned this way.** `ActsExamples::Geant4Simulation`
  forwards `simGeoContext` to the user actions, but the G4 geometry itself is
  built once and is not context aware, so per-event transforms never reach G4
  transport. Use Fatras.
- **Navigation is built on the nominal geometry.** Layer arrays, volume
  boundaries and surface binning are constructed once, so misalignments have to
  stay small enough not to break navigation.

## Where to look

- `Examples/Scripts/Python/misaligned_simulation.py` - a telescope with one
  shifted layer, showing the residual bias that appears when only simulation
  sees the shift.
- `Examples/Scripts/Python/millepede_alignment.py` - the same setup feeding a
  Millepede alignment that fits the shift back out.
- `Python/Examples/tests/test_alignmentdecorator.py` - the corresponding tests.
