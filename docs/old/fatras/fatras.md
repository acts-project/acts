# Fatras

Fatras is an ACTS implementation for fast track simulation which was primarily built to easily obtain hits from a tracking geometry.
Note that it is not meant for detailed physical simulations or studies.

Fatras is using the {class}`Acts::Propagator` as a backend for the track extrapolation. The simulation is just an **Actor** in the terminology of the ACTS propagation.

:::{tip}
A dedicated description of the propagation can be found [here](propagation_impl).
:::

Hits on sensitive surfaces are recorded and tagged with the particle passing through the surface.
Later, this hit information can be used by a digitization algorithm to mimic a detector response.

ACTS Fatras is fully capable of handling secondary particles even though the current set of interactions is not utilizing this to a big extent.
The simulation will propagate each particle until it reaches the end of the detector or until a specified path length is reached.
Afterwards, the first secondary particle that might have been generated in the process will be propagated.
This continues until all particles are transported.

## Barcode identifiers

Fatras labels simulated particles and hits with `ActsFatras::Barcode`.
A barcode is the event-local particle identifier used by Fatras and by the ACTS examples framework for truth matching.
It is not a geometry identifier; it answers which simulated particle produced a hit, particle state, or truth-matching entry.

The barcode stores five integer components:

| Component | Meaning |
| --------- | ------- |
| primary vertex | The primary interaction vertex. Ordinary simulated particles use a non-zero primary vertex. |
| secondary vertex | A secondary vertex below the primary vertex. Zero means the particle comes directly from the primary vertex. |
| particle | The particle number within the selected primary and secondary vertex. |
| generation | The descendant generation. Particles produced at the vertex use generation zero. |
| sub-particle | The particle number within a non-zero generation. Particles produced at the vertex use sub-particle zero. |

The default-constructed barcode has all components set to zero and represents an invalid, missing, or unknown particle identifier.
This value is useful as a sentinel, but it should not be used for ordinary simulated particles.

For example, a particle from primary vertex `2` with particle number `14` is encoded as:

```text
vp=2|vs=0|p=14|g=0|sp=0
```

If a Fatras interaction creates two descendant particles from this particle, the new particles keep the same vertex and particle number.
The generation is increased and the sub-particle number distinguishes the two descendants:

```text
vp=2|vs=0|p=14|g=1|sp=0
vp=2|vs=0|p=14|g=1|sp=1
```

The preserved vertex and particle components make it cheap to recover the initial simulated particle for truth matching, while the generation and sub-particle components distinguish particles created during the simulation.
Interactions that create descendants can use `ActsFatras::Barcode::makeDescendant` to increase the generation and set the sub-particle number.

Two helper projections are commonly useful when grouping truth information:

- `ActsFatras::Barcode::vertexId` drops the particle and sub-particle components so objects can be grouped by production vertex.
- `ActsFatras::Barcode::withoutSubparticle` drops only the sub-particle component so generated descendants can be grouped by original particle and generation.

Since barcodes are created locally by the code that produces particles, there is no global allocation service that stores the full decay tree.
Independent interactions can therefore create overlapping descendant identifiers if they start from the same particle and generation.
When the full set of particles is available, sub-particle numbers within a generation can be renumbered to make the identifiers unique.

## Supported interactions

Fatras implements a few interactions
 - Bethe-Bloch energy loss for charged particles (see `ActsFatras::BetheBloch`)
 - Bethe-Heitler energy loss for electrons (and positrons) (see `ActsFatras::BetheHeitler`)
 - Photon conversion for pair production (see `ActsFatras::PhotonConversion`)
 - Scattering for charged particles (see `ActsFatras::GenericScattering`)

These interactions are meant to be physical accurate within their boundaries but are far from sufficient for a full simulation.

## Use-cases

The primary use-case for ACTS Fatras is early algorithm validation.
It is fast enough to simulate thousands of small to medium events (about 100 tracks) on consumer hardware in a few seconds.
This provides a quick turnaround time for development.

Fatras does not replace a full detector simulation like Geant4 as only a few interactions are implemented and the geometry is highly simplified.
