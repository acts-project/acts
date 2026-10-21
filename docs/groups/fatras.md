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
