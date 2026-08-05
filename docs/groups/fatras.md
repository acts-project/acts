@defgroup fatras Fatras
@brief Fast track simulation on ACTS tracking geometries.

# Fatras

Fatras is the ACTS fast track simulation package. It uses the
@ref Acts::Propagator and a reconstruction-oriented @ref Acts::TrackingGeometry
to transport particles through the detector and to produce simulated hits on
sensitive surfaces. It is intended for fast algorithm validation and detector
studies that do not require a full, detailed detector simulation.

Fatras uses parametrized material interactions, including multiple scattering,
Bethe-Bloch energy loss, Bethe-Heitler energy loss, and photon conversion. Since
it operates on tracking surfaces rather than a full volumetric detector model,
it is faster than a detailed Geant4 simulation but also less complete.

## Synthetic space point events {#fatras_synthetic_events}

@ref ActsFatras::Synthetic is a second, much smaller fast simulation that shares
none of the machinery above. It fills an @ref Acts::SpacePointContainer directly
-- no @ref Acts::TrackingGeometry, no @ref Acts::Propagator, no input data -- and
makes an ATLAS-like pile-up of 200 in tens of milliseconds. It answers one
question cheaply: what does a seeder see on a busy event? Space points carry the
layer and the particle they came from, so that needs no truth matching.

It is deliberately coarse: a tool for throughput and combinatorics studies rather
than for physics performance.

### Geometry

A @ref ActsFatras::Synthetic::DetectorLayout is cylinders at fixed `r` and discs
at fixed `z`, as plain structs, so a helix crosses every surface in closed form.
Nothing is resolved in azimuth and a crossing leaves at most one hit. The
restriction is on the *shape* alone -- every cylinder has its own half-length and
every disc its own rings -- so an endcap can be described as the staggered discs
and rings it is.

- @ref ActsFatras::Synthetic::DetectorLayoutBuilder builds one from scratch;
  descriptions of the ATLAS ITk, ODD and ACTS Generic pixels ship with it.
- @ref ActsFatras::Synthetic::PassiveSurfaceDescription adds the supports and
  services a layer is carried on, which produce secondaries like any material.
- @ref ActsFatras::Synthetic::SurfaceMaterial bands a surface along the
  coordinate it extends in, so a ring is told from the support beside it. The
  shipped layouts carry what the detector's own material map reports.

### Event content

@ref ActsFatras::Synthetic::EventGenerator draws the primaries -- a rapidity
plateau with a Fermi edge, Gaussian in `z0` and `d0`, with a Hagedorn spectrum
-- then walks each track through the layout and makes what it makes.

What a crossing does, each switchable in @ref ActsFatras::Synthetic::EventConfig:

- **Secondaries**, at the material the surface carries times the path length
  through it and nothing fitted on top. Electrons count per `X0` and nuclear
  products per `L0`, in three channels: a radial electron, a forward cascade
  product, and an isotropic evaporation product.
- **Stubs**, where a daughter is too soft to leave the surface that made it.
- **Multiple scattering and energy loss**, which displace a hit and leave the
  trajectory alone. This is what gives a seed the spread a seeder cuts on. A
  track that cannot pay for the surface in front of it ranges out.
- **Module overlaps**, adjacent modules alternating in depth so that a track
  through their common edge is measured twice.

A track is followed past its outermost point, so a soft one curls back through
the layers it has already crossed, bounded by
@ref ActsFatras::Synthetic::PropagationConfig::maxTurns and by the escape bounds
of the enclosing tracker. Neutral long-lived particles decaying near the beam
line are the only secondaries produced away from a surface.

### Tuned configurations

The detector enters the beam spot, the resolution and both yields, so a
configuration belongs to a layout. @ref ActsFatras::Synthetic::EventConfig ships
@ref ActsFatras::Synthetic::EventConfig::itkPixelTtbarPu200, fitted against a
GNN4ITk Athena dump, and
@ref ActsFatras::Synthetic::EventConfig::openDataDetectorTtbarPu200 against
ColliderML. Positions and material are read off the detector description and the
overlaps and secondary kinematics measured on the reference; only the yields, the
spectrum and the beam spot are fitted.

Each is fitted on one half of its sample and checked on the other. On the
held-out half, per event and normalised to the reference:

| | ITk | ODD |
| --- | --- | --- |
| space points | 0.99 | 0.99 |
| &nbsp;&nbsp;primary, inside the generated acceptance | 0.97 | 1.02 |
| &nbsp;&nbsp;primary, outside it | 1.03 | 1.16 |
| &nbsp;&nbsp;non-primary | 1.00 | 0.95 |
| primaries / event | 0.98 | 1.00 |
| mean primary hits | 0.99 | 1.02 |
| mean secondary hits | 0.97 | 1.10 |

Known to be off: the ODD's secondary momentum runs a fifth high and its `|d0|`
reach short, its forward production is short beyond `|z| = 900`, and the spectrum
is over-produced below the 100 MeV the references stop recording at.

Counting secondaries takes care -- a full simulation records one only above a
truth-link threshold, and two thirds of the real ones fall below it. Compare
non-primary space points rather than particle counts.

## Barcode identifiers {#fatras_barcode_identifiers}

Fatras labels simulated particles and hits with @ref ActsFatras::Barcode. A
barcode is the event-local particle identifier used by Fatras and by the ACTS
examples framework for truth matching. It is not a geometry identifier; it
answers "which simulated particle produced this object?".

The barcode stores five integer components:

| Component | Meaning |
| --------- | ------- |
| primary vertex | The primary interaction vertex. Ordinary simulated particles use a non-zero primary vertex. |
| secondary vertex | A secondary vertex below the primary vertex. Zero means the particle comes directly from the primary vertex. |
| particle | The particle number within the selected primary and secondary vertex. |
| generation | The descendant generation. Particles produced at the vertex use generation zero. |
| sub-particle | The particle number within a non-zero generation. Particles produced at the vertex use sub-particle zero. |

The default-constructed barcode has all components set to zero and represents an
invalid, missing, or unknown particle identifier. This value is useful as a
sentinel, but it should not be used for ordinary simulated particles.

For example, a particle from primary vertex `2` with particle number `14` is
encoded as:

```text
vp=2|vs=0|p=14|g=0|sp=0
```

If a Fatras interaction creates two descendant particles from this particle, the
new particles keep the same vertex and particle number. The generation is
increased and the sub-particle number distinguishes the two descendants:

```text
vp=2|vs=0|p=14|g=1|sp=0
vp=2|vs=0|p=14|g=1|sp=1
```

The preserved vertex and particle components make it cheap to recover the
initial simulated particle for truth matching, while the generation and
sub-particle components distinguish particles created during the simulation.

## Creating and reducing barcodes

@ref ActsFatras::Barcode is immutable from the caller's point of view: modifier
methods return a new barcode with one component changed. Typical construction
therefore chains the `with...` methods:

```cpp
auto particleId = ActsFatras::Barcode()
                      .withVertexPrimary(1)
                      .withVertexSecondary(0)
                      .withParticle(42);
```

Interactions that create descendants can call
@ref ActsFatras::Barcode::makeDescendant to increase the generation and set the
sub-particle number:

```cpp
auto electronId = photonId.makeDescendant(0);
auto positronId = photonId.makeDescendant(1);
```

Two helper projections are commonly useful when grouping truth information:

- @ref ActsFatras::Barcode::vertexId drops the particle and sub-particle
  components so objects can be grouped by production vertex.
- @ref ActsFatras::Barcode::withoutSubparticle drops only the sub-particle
  component so generated descendants can be grouped by original particle and
  generation.

Since barcodes are created locally by the code that produces particles, there is
no global allocation service that stores the full decay tree. Independent
interactions can therefore create overlapping descendant identifiers if they
start from the same particle and generation. When the full set of particles is
available, sub-particle numbers within a generation can be renumbered to make
the identifiers unique.
