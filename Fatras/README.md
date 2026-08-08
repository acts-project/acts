# FAst TRAck Simulation package

This packages provides tools to run a fast track simulation on top of the core
library. The fast track simulation (Fatras) uses the actor plug-in mechanism of
the `Acts::Propagator` and its predictive navigation to simulate particle
trajectories through the tracking detector.

As a fast track simulation, it uses the surface-based reconstruction geometry
description, i.e. `Acts::TrackingGeometry`, as a simulation geometry instead of
a detailed volumetric description. Interactions and material effects are
simulated using parametrized models:

-   Multiple Coulomb scattering is simulated by Gaussian (mixture)
    approximations.
-   Ionisation loss is simulated using the Bethe-Bloch formalism.
-   Radiation loss is simulated using the Bethe-Heitler formalism.

## Synthetic space point events

The `ActsFatras::Synthetic` subcomponent is the exception to everything above: it
uses neither the propagator nor `Acts::TrackingGeometry`. Its geometry is a list
of cylinders around the beam axis and discs normal to it as plain structs, i.e. a
collider-style barrel and endcap detector, so a helix crosses every surface in
closed form, and its output is an `Acts::SpacePointContainer` rather than hits on
surfaces. That makes it fast enough to generate an ATLAS-like pile-up of 200 in
tens of milliseconds, which is what a seeding benchmark needs and what the
propagator-based simulation above cannot provide.

It is deliberately coarse in its physics and is meant for throughput and
combinatorics studies rather than physics performance. Its *geometry* is not:
every disc holds its own rings and every barrel cylinder its own half-length, so
the ITk pixel endcap is described as the seventy-five staggered discs and
ninety-five rings it is.

An event configuration is fitted against a full simulation of one detector *and*
one layout of it, so a layout should be used with the configuration of the same
name. `docs/groups/fatras.md` documents what such a fit reproduces and what it
cannot.
