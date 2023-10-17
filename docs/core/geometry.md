(geometry_impl)=
# Geometry module

The ACTS geometry model is strongly based on the ATLAS Tracking geometry. Its
core is built on a surface-based description that make up all geometry objects
of higher complexity. This design has been chosen as the surface objects can be
used together with the track propagation module and thus all geometry objects
become natively integrated into the tracking software.

```{note}
There is an ongoing rewrite of the geometry and navigation modules where logical layers will be modelled as volumes, see [](exp_geometry_impl).

```

## GeometryObject base class

All geometry objects in Acts inherit from a virtual {class}`Acts::GeometryObject` base class.

:::{doxygenclass} Acts::GeometryObject
:::


```{note}
The `binningPosition(const GeometryContext& gctx, BinningValue& bVal)` method allows to define 
a customized position of this object when being subject to ordering in some specific binnings, e.g. in a `BinnedArray` or a `SurfaceArray`.
```

This class ensures that a unique {class}`Acts::GeometryIdentifier` is assigned to every geometry
object. 

## Geometry identifier

The {class}`Acts::GeometryIdentifier` is mainly used for fast identification of the type of
the geometry object (as most of them are either extensions or containers of the
{class}`Acts::Surface` objects) and for the identification of the geometry surfaces after
building, e.g. for the uploading/assigning of material to the surface after
creation. The {class}`Acts::GeometryIdentifier` uses a simple masking procedure for applying an
identification schema.

While it is used in Acts-internal applications such as material mapping, it is not employed for
`EventData` and `Geometry` identification in an experiment setup. Instead, one should define and use the
`Identifier` class in the latter case.

:::{doxygenclass} Acts::GeometryIdentifier
---
members: kVolumeMask,kBoundaryMask,kLayerMask,kApproachMask,kSensitiveMask,kExtraMask
---
:::


## Surface classes

All classes which represent a thin surface in ACTS inherit from
the common virtual base class {class}`Acts::Surface`, which defines
the public interface of all surfaces. While the different concrete
surface classes are defined by their respective native local
coordinate system, the shapes on these surfaces are defined by classes
that inherit from {class}`Acts::SurfaceBounds`, which every surface must provide.
In case of boundless surfaces, a special {class}`Acts::InfiniteBounds` class is
available.

Each {class}`Acts::Surface` instance reports its type from {func}`Acts::Surface::type()`:

:::{doxygenenum} Acts::Surface::SurfaceType
:::


| Surface Type          | Local Coordinates | Bound Types available                                                                       |
|:----------------------|-------------------|:--------------------------------------------------------------------------------------------|
| {class}`Acts::ConeSurface`         | $[r\phi, z]$         | {class}`Acts::ConeBounds`                                                                                |
| {class}`Acts::CylinderSurface`     | $[r, \phi]$         | {class}`Acts::CylinderBounds`                                                                            |
| {class}`Acts::DiscSurface`         | $[r, \phi]$         | {class}`Acts::RadialBounds`, {class}`Acts::DiscTrapezoidalBounds`                                                     |
| {class}`Acts::PlaneSurface`        | $[x, y]$         | {class}`Acts::RectangleBounds`, {class}`Acts::TrapezoidalBounds`, <br>{class}`Acts::TriangleBounds`,{class}`Acts::InfiniteBounds`, <br> {class}`Acts::EllipseBounds` |
| {class}`Acts::PerigeeSurface`,<br> {class}`Acts::StrawSurface` | $[d, z]$ | {class}`Acts::CylinderBounds`                                                                            |
| {class}`Acts::LineSurface` | $[d_0, z_0]$ | {class}`Acts::LineBounds` |

```{tip}
In an ideal setup, the coordinate systems also define the readout
measurement directions. In such a case, a track prediction from the
propagation will already be in the correct frame of the measurement and
residual or compatibility checks will not need additional coordinate
transformations.
```

### Plane surface

![PlaneBounds](/figures/geometry/PlaneBounds.png)

:::{doxygenclass} Acts::PlaneSurface
---
members: globalToLocal,localToGlobal,intersect,normal
---
:::

### Disc surface

![DiscBounds](/figures/geometry/DiscBounds.png)

:::{doxygenclass} Acts::DiscSurface
---
members: globalToLocal,localToGlobal,intersect,normal
---
:::
 
### Cylinder surface

![CylinderBounds](/figures/geometry/CylinderBounds.png)

:::{doxygenclass} Acts::CylinderSurface
---
members: globalToLocal,localToGlobal,intersect,normal
---
:::

### Cone surface

:::{doxygenclass} Acts::ConeSurface
---
members: globalToLocal,localToGlobal,intersect,normal
---
:::

### Line surface

{class}`Acts::LineSurface` is a special kind of surface that depends on a reference
direction, typically the unit momentum direction $\vec d$ of a particle. A point in
space is considered *on surface* if and only if it coincides with the point of
closest approach between the direction vector $\vec d$ and the line direction
vector $\vec z$. As such, the function {func}`Acts::LineSurface::globalToLocal`
can fail, if the argument position and direction do not fulfill this criterion.
It is pure-virtual, meaning that it can not be instantiated on its own.

:::{doxygenclass} Acts::LineSurface
---
members: globalToLocal,localToGlobal,intersect,normal
---
:::

#### Straw surface

:::{doxygenclass} Acts::StrawSurface
---
members: false
---
:::

#### Perigee surface

:::{doxygenclass} Acts::PerigeeSurface
---
members: false
---
:::


## Layer classes


The {class}`Acts::Layer` class is an extension of the {class}`Acts::Surface` class that allows the
definition of sub surfaces (sensitive surfaces for modules, or extra material
surfaces).

The layer can simply correspond to a 'virtual' surface in the detector
description or represent a more complex object that may contain:

* a representing surface, which is accessible via a {func}`Acts::Layer::surfaceRepresentation`
* an array of contained surfaces, accessible via {func}`Acts::Layer::surfaceArray` method
* approach surfaces (i.e. boundary surface of the volume occupied by the layer)
* surface material description on any of the confined surfaces

The following illustration shows an $xy$ view of a cylinder layer with planar
detection modules:

![CylinderLayer](/figures/geometry/CylinderLayer.png)

Modules can be sorted onto layer using all supported binning methods described
through the {class}`Acts::SurfaceArray` class. The binning can be adjusted to fit as well as
possible.

![DiscLayerEB](/figures/geometry/DiscLayerEB.png)

The unoccupied space in a volume that contains a layer array is filled with
objects of type {class}`Acts::NavigationLayer`, which allows that in a fully static geometry
setup, every single point in a volume can be associated with a layer. Layer
objects are confined together in a special {type}`Acts::LayerArray` class and can be
contained by a {class}`Acts::TrackingVolume`.

![LayerArray](/figures/geometry/LayerArray.png)



## Volume classes

The {class}`Acts::Volume` class is a container of
{type}`Acts::BoundarySurface` objects, where each
{type}`Acts::BoundarySurface` is an extension of the {class}`Acts::Surface`
class with additional information about the attached volumes. The normal vector
of the surface defines an *inside* (opposite w.r.t. the normal vector) and an
*outside* (along w.r.t. the normal vector) direction. Either a single volume or
an array of volumes can be attached to a volume.

The simples volume class is just a collection of surfaces, where the
{class}`Acts::TrackingVolume` describes a volume that can contain:

* an array of contained layers
* an array of contained volumes (as a container volume)
* an array of contained volumes (as *floating* objects)
* a volume based material description

The shape of the volume is defined by {class}`Acts::VolumeBounds` classes that create the
corresponding bounding surfaces and register the attachment to the volume itself
at creation.

![VolumeBounds](/figures/geometry/VolumeBounds.png)
![CylinderVolumeBounds](/figures/geometry/CylinderVolumeBounds.png)

## Detector material description

Two types of material description exist, one for a surface based material, one
for a volume based material. They will be dealt with differently in the
extrapolation.

The basic information for any material is:

* the radiation length X0 the nuclear interaction length L0 the atomic weight A
* the atomic charge Z the density of the material

This information is confined together in the {class}`Acts::Material` class.

```{note}
In track reconstruction, only an effective material description is needed, i.e.
non-physical values in regards of the atomic number, the elementary charge or
even the density are allowed, as long as the effective relative radiation
length and $A/Z \times \rho$ ratio can be retrieved.  This enables the compactification
of the material description, as the element composition record does not have to
be kept.
```

Surface based material extends this material information by a representative
thickness; the corresponding object is called {class}`Acts::MaterialSlab`. The
thickness hereby can be arbitrarily chosen in order to regulate the material
budget, it does not have to represent the actual thickness of a detector
element. To attach it to a surface, a dedicated {class}`Acts::SurfaceMaterial`
class (or it's extensions) is used, which allows to also describe binned
material.

Possible extensions are:

 * {class}`Acts::HomogeneousSurfaceMaterial`, homogeneous material description on a surface
 * {class}`Acts::BinnedSurfaceMaterial`, an arbitrarily binned material description with a
    corresponding {class}`Acts::BinUtility` object

In addition, a dedicated extension exists to allow configuration of the material
mapping process, that is in further described below.

 * {class}`Acts::ProtoSurfaceMaterial`, only binning description (without material) to be
   used in the material mapping process

## Geometry building

The geometry building procedure follows the ATLAS tracking geometry philosophy of
a static frame of *glued* volumes, that lead the navigation flow through the
geometry,

### Attaching a 3D detector geometry

Usually, a 3D detector model geometry exists, which is either native to the full
detector simulation (Geant4) or is translated into it. This model, however, is
in general too detailed for track reconstruction: navigating through the
detailed detector geometry is generally costly and one can profit greatly from a simplification mechanism.

For most part of the track reconstruction, only a surface based description of
the detector is needed, in order to allow (surface based) material integration
and parametrization/prediction of trajectories on detection surfaces. It is thus
necessary that the detection surfaces are described to full detail in the
reconstruction geometry (called {class}`Acts::TrackingGeometry`). This is guaranteed by a
proxy mechanism that connects the detection elements (conveniently called
{class}`Acts::DetectorElement`) to {class}`Acts::Surface` object in the reconstruction:

![DetectorElement](/figures/geometry/DetectorElement.png)

#### Existing plugins for 3D geometry libraries

Very simple helper methods for 3D libraries exist, they are certainly not
optimised, but used for templating:

* {class}`Acts::TGeoDetectorElement` connects a TGeo volume to a {class}`Acts::Surface`
* {class}`Acts::DD4HepDetectorElement` connects a DD4hep volume (based on TGeo) to a {class}`Acts::Surface`
* {class}`Acts::Geant4DetectorElement` connects a Geant4 volume to a {class}`Acts::Surface`

Further extensions exist in dedicated experiment contexts, such as e.g. a `GeoModel`
binding for the ATLAS experiment.

```{note}
While `DD4hep` offers a descriptive language with a dedicated extension mechanism
that can be used by Acts to interpret the underlying geometry hierarchy and and structure,
there is no such guarantee when having the already as built `TGeo` geometry in hand.
Therefore a dedicated Acts configuration file based on `json` can be provided that allows
to specify parsing restrictions for sub detectors. 
```


### Layer building

{class}`Acts::Surface` objects that are to be grouped on a layer should be put into a
{class}`Acts::SurfaceArray` and provided to the layer. Certain helper tools exist to ease the
translation and create appropriate binning structure: The {class}`Acts::SurfaceArrayCreator`
can create cylindrical, disc-like & planar layers, where the dimensions of the
layer are determined by parsing the provided surfaces. Additionally, an envelope
covering the surfaces can be chosen.

```{note}
There exist standard layer builders that are designed to build cylindrical, disk like 
and planar layers and perform the ordering of the surfaces onto those layers. These
builders are called from the top level translation entry points from either `TGeo` 
or `DD4hep`.
```


### Volume building, packing, and gluing

The philosophy of the {class}`Acts::TrackingGeometry` is a fully connective geometry setup,
i.e. {class}`Acts::TrackingVolume` objects are either pure containers for other contained
{class}`Acts::TrackingVolume` instances (where the contained volumes fully fill the space of
the container volume), or are fully attached via the boundary surface mechanism.
The boundary surfaces then act as portals from one {class}`Acts::TrackingVolume` into the
next one along the trajectory.

The process to create a fully connected tracking geometry is called glueing.
Wherever possible, common boundary surfaces are *shared*, where this is not
possible, they are *attached*.

![GlueingBC](/figures/geometry/GlueingBC.png)
![GlueingABC](/figures/geometry/GlueingABC.png)
![NavigationABC](/figures/geometry/NavigationABC.png)

For cylindrical detector setups, a dedicated {class}`Acts::CylinderVolumeBuilder` is
provided, which performs a variety of volume building, packing and gluing.

```{note}
For most cylindrical detectors, there exist automated glueing and geometry building 
modules that take care of the glueing process.
```

## TrackingGeometry building using a KDTree and a Proto Description

For cylindrical detectors there exist a generic tracking geometry building module,
based on KDTree and a proto description.

This building procedure uses a {class}`Acts::ProtoDetector` description which provides a 
high level description of layers and container volumes, together with some 
binning and ordering information.
This proto description is then used to assign surfaces that are provided to the 
{class}`Acts::KDTreeTrackingGeometryBuilder` using an internal query to the KD-tree structure.
