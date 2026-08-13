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

A detector is written down as a @ref ActsFatras::Synthetic::DetectorDescription
and expanded into that layout by @ref ActsFatras::Synthetic::makeLayout. The
description is the form everything else works on, and it is a set of named
@ref ActsFatras::Synthetic::SubsystemDescription -- the pixels, the strips, one of
a detector's systems -- each fanning out into barrels and endcaps. A name is what
makes a subsystem selectable and what material is keyed onto, so one detector can
be built whole or a system at a time.

- @ref ActsFatras::Synthetic::CylinderDescription and
  @ref ActsFatras::Synthetic::DiscDescription are a barrel's and an endcap's
  layers. Each carries its own material, module structure and overlap numbers,
  and answers to a layer index of its own rather than to its position in a list.
- @ref ActsFatras::Synthetic::EndcapPlacement lets an endcap say which sides of
  the interaction point it is built on, positions being quoted as absolute `z`,
  so a symmetric detector is written once and an asymmetric one can still be
  spelled out.
- @ref ActsFatras::Synthetic::PassiveSurfaceDescription adds the supports and
  services a layer is carried on, which produce secondaries like any material.
  The beam pipe is one of these and belongs to the detector rather than to any
  subsystem, being the only material in front of the innermost layer.
- @ref ActsFatras::Synthetic::SurfaceMaterial bands a surface along the
  coordinate it extends in, so a ring is told from the support beside it. The
  bands are the gaps between their edges and follow where the material actually
  changes, not the rings; the shipped descriptions carry what the detector's own
  material map reports.
- @ref ActsFatras::Synthetic::DetectorLayoutBuilder assembles a layout surface by
  surface where a description would be in the way. A description of the ACTS
  Generic pixels ships in C++, and is the one detector that does.

A detector is held whole and built in parts:
@ref ActsFatras::Synthetic::selectSubsystems narrows it to the systems named,
keeping the beam pipe and the containment of the whole tracker whichever they
are, and @ref ActsFatras::Synthetic::merge puts descriptions back together, so a
hand-written subsystem can be added to a shipped detector. A space point still
says which system it came from afterwards, its layer carrying an index into
@ref ActsFatras::Synthetic::DetectorLayout::subsystems.

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

### Configuration

The detector enters the beam spot, the resolution and both yields, so a
configuration belongs to a layout: positions and material are read off the
detector description and the overlaps and secondary kinematics measured on a
full simulation of it, leaving the yields, the spectrum and the beam spot to be
fitted. The defaults of @ref ActsFatras::Synthetic::EventConfig belong to no
detector in particular.

Counting secondaries to fit against takes care -- a full simulation records one
only above a truth-link threshold, and two thirds of the real ones fall below
it. Compare non-primary space points rather than particle counts.

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
