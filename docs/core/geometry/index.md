(geometry_impl)=
# Geometry module

The ACTS geometry model is strongly based on the ATLAS Tracking geometry. Its
core is built on a surface-based description that make up all geometry objects
of higher complexity. This design has been chosen as the surface objects can be
used together with the track propagation module and thus all geometry objects
become natively integrated into the tracking software.

```{note}
There is an ongoing rewrite of the geometry and navigation modules where
logical layers will be modelled as volumes, see [](layerless_geometry).

```

:::{toctree}
:maxdepth: 1
geometry_id
material
surfaces
legacy/legacy
layerless/layerless
:::
