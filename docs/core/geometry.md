(geometry_impl)=
# Geometry module

The ACTS geometry model is strongly based on the ATLAS Tracking geometry. Its
core is built on a surface-based description that make up all geometry objects
of higher complexity. This design has been chosen as the surface objects can be
used together with the track propagation module and thus all geometry objects
become natively integrated into the tracking software.

## GeometryObject base class and GeometryIdentifier

All geometry objects in Acts inherit from a virtual `GeometryObject` base class

```cpp
class GeometryObject {
 public:
  /// Defaulted construrctor
  GeometryObject() = default;

  /// Defaulted copy constructor
  GeometryObject(const GeometryObject&) = default;

  /// Constructor from a value
  ///
  /// @param geometryId the geometry identifier of the object
  GeometryObject(const GeometryIdentifier& geometryId)
      : m_geometryId(geometryId) {}

  /// Assignment operator
  ///
  /// @param geometryId the source geometryId
  GeometryObject& operator=(const GeometryObject& geometryId) {
    if (&geometryId != this) {
      m_geometryId = geometryId.m_geometryId;
    }
    return *this;
  }

  /// @return the geometry id by reference
  const GeometryIdentifier& geometryId() const;

  /// Force a binning position method
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the value in which you want to bin
  ///
  /// @return vector 3D used for the binning schema
  virtual Vector3 binningPosition(const GeometryContext& gctx,
                                  BinningValue bValue) const = 0;

  /// Implement the binningValue
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the dobule in which you want to bin
  ///
  /// @return float to be used for the binning schema
  virtual double binningPositionValue(const GeometryContext& gctx,
                                      BinningValue bValue) const;

  /// Set the value
  ///
  /// @param geometryId the geometry identifier to be assigned
  void assignGeometryId(const GeometryIdentifier& geometryId);

 protected:
  GeometryIdentifier m_geometryId;
};
```

This class ensures that a unique `GeometryIdentifier` is assigned to every geometry
object. The `GeometryIdentifier` is mainly used for fast identification of the type of
the geometry object (as most are either extensions or containers of the
`Surface` objects) and for the identification of the geometry surfaces after
building, e.g. for the uploading/assigning of material to the surface after
creation. The `GeometryIdentifier` uses a simple masking procedure for applying an
identification schema.

It is used for Acts internal applications, such as material mapping, but not for
`EventData` and `Geometry` identification in an experiment setup, for this the
`Identifier` class is to be used and/or defined.

```cpp
typedef uint64_t geo_id_value;

namespace Acts {

/// @class GeometryIdentifier
///
///  Identifier for Geometry nodes - packing the
///  - (Sensitive) Surfaces    - uses counting through sensitive surfaces
///  - (Approach)  Surfaces    - uses counting approach surfaces
///  - (Layer)     Surfaces    - uses counting confined layers
///  - (Boundary)  Surfaces    - uses counting through boundary surfaces
///  - Volumes                 - uses counting given by TrackingGeometry
class GeometryIdentifier
{
public:
  const static geo_id_value volume_mask    = 0xff00000000000000;
  const static geo_id_value boundary_mask  = 0x00ff000000000000;
  const static geo_id_value layer_mask     = 0x0000ff0000000000;
  const static geo_id_value approach_mask  = 0x000000f000000000;
  const static geo_id_value sensitive_mask = 0x0000000ffff00000;
  const static geo_id_value channel_mask   = 0x00000000000fffff;
  ...
};
```

```{note}
The `binningPosition(const GeometryContext& gctx, BinningValue& bVal)` method allows to define 
a customized position of this object when being subject to ordering in some specific binnings, e.g. in a `BinnedArray` or a `SurfaceArray`.
```

## Surface classes

The `Surface` class builds the core class of all geometry objects and can be
used natively with the propagation and extrapolation modules. The common
`Surface` virtual base defines the public interface of all surfaces. The
different concrete `Surface` classes are defined by their respective native
local coordinate system, while different shapes on surfaces are defined by
`SurfaceBounds` classes which every surface must provide. In case of boundless
surfaces, a special `InfiniteBounds` class is available.

| Surface Type          | Local Coordinates | Bound Types available                                                                       |
|:----------------------|-------------------|:--------------------------------------------------------------------------------------------|
| `ConeSurface`         | [rphi, z]         | `ConeBounds`                                                                                |
| `CylinderSurface`     | [r, phi]          | `CylinderBounds`                                                                            |
| `DiscSurface`         | [r, phi]          | `RadialBounds`, `DiscTrapezoidalBounds`                                                     |
| `PlaneSurface`        | [x, y]            | `RectangleBounds`, `TrapezoidalBounds`, `TriangleBounds`,<br> `InfiniteBounds`, `EllipseBounds` |
| `PerigeeSurface`, `StrawSurface` | [d, z] | `CylinderBounds`                                                                            |

![CylinderBounds](/figures/geometry/CylinderBounds.png)
![DiscBounds](/figures/geometry/DiscBounds.png)
![PlaneBounds](/figures/geometry/PlaneBounds.png)

```{tip}
The coordinate systems define - in an ideal setup - also the readout measurement directions. In such a case, a track prediction from the propagation, e.g. will already be in the correct frame of the measurement and residual or compatibility checks will not need additional coordinate transformations.
```

## Layer classes

The `Layer` class is an extension of the `Surface` class that allows the
definition of sub surfaces (sensitive surfaces for modules, or extra material
surfaces).

The Layer can simply correspond to a 'virtual' surface in the detector
description or represent a more complex object that may contain:

* a representing surface, which is accessible via a `representingSurface()`
* method an array of contained surfaces, accessible via `surfaceArray()` method
* approach surfaces (i.e. boundary surface of the volume occupied by the layer)
* surface material description on any of the confined surfaces

The following illustration shows an x-y view of a cylinder layer with planar
detection modules:

![CylinderLayer](/figures/geometry/CylinderLayer.png)

Modules can be sorted onto layer using all supported binning methods described
through the `SurfaceArray` class, the binning can be adjusted to fit as good as
possible.

![DiscLayerEB](/figures/geometry/DiscLayerEB.png)

The un-occupied space in a volume which contains a layer array is filled with
objects of type `NavigationLayer`, which allows that in a fully static geometry
setup, every single point in a volume can be associated with a layer. Layer
objects are confined together in a special `LayerArray` class and can be
contained by a `TrackingVolume`.

![LayerArray](/figures/geometry/LayerArray.png)


## Volume classes

The `Volume` class is a container of `BoundarySurface` objects, where each
`BoundarySurface` is an extension of the `Surface` class with additional
information about the attached Volumes. The normal vector of the surface
defines an *inside* (opposite w.r.t. the normal vector) and an *outside* (along
w.r.t. the normal vector) direction. Either a single volume or an array of
volumes can be attached to a volume.

The simples volume class is just a collection of surfaces, where the
`TrackingVolume` describes a volume that can contain:

* an array of contained layers
* an array of contained volumes (as a container volume)
* an array of contained volumes (as *floating* objects)
* a volume based material description

The shape of the volume is defined by `VolumeBounds` classes that create the
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

This information is confined together in the `Material` class.

```{note}
In track reconstruction, only an effective material description is needed, i.e. non-physical values in regards of the atomic number, the elementary charge or even the density are allowed, as long as the effective relative radiation length and `A/Z x rho` ratio can be retrieved.
This enables compatictification of the material description, as the element composition record does not have to be kept.
```

Surface based material extends this material information by representative
thickness, the corresponding object is called `MaterialSlab`. The
thickness hereby can be arbitrarily chosen in order to regulate the material
budget, it does not have to represent the actual thickness of a detector
element. To attach it to a surface, a dedicated `SurfaceMaterial` class (or it's
extensions) is used, which allows to also describe binned material.

Possible extensions are:

 * `HomogeneousSurfaceMaterial`, homogeneous material description on a surface
 * `BinnedSurfaceMaterial`, an arbitrarily binned material description with a
    corresponding `BinUtility` object

In addition, a dedicated extension exists to allow configuration of the material
mapping process, that is in further described below.

 * `ProtoSurfaceMaterial`, only binning description (without material) to be
   used in the material mapping process

## Geometry building

The geometry building procedure follows the ATLAS TrackingGeometry philosophy of
a static frame of *glued* volumes, that lead the navigation flow through the
geometry,

### Attaching a 3D detector geometry

Usually, a 3D detector model geometry exists, which is either native to the full
detector simulation (Geant4) or is translated into it. This model, however, is
in general too detailed for track reconstruction: navigating through the
detailed detector geometry is generally costly and one can profit greatly from a simplication mechanism.

For most part of the track reconstruction, only a surface based description of
the detector is needed, in order to allow (surface based) material integration
and parametrization/prediction of trajectories on detection surfaces. It is thus
necessary that the detection surfaces are described to full detail in the
reconstruction geometry (called `TrackingGeometry`). This is guaranteed by a
proxy mechanism that connects the detection elements (conveniently called
`DetectorElement`) to `Surface` object in the reconstruction:

![DetectorElement](/figures/geometry/DetectorElement.png)

#### Existing plugins for 3D geometry libraries

Very simple helper methods for 3D libraries exist, they are certainly not
optimised, but used for templating:

* `TGeoDetectorElement` connects a TGeo volume to a `Surface`
* `DD4HepDetectorElement` connects a DD4hep volume (based on TGeo) to a `Surface`
* `Geant4DetectorElement` connects a Geant4 volume to a `Surface`

Further exensions exist in dedicated experiment contexts, such as e.g. a `GeoModel`
binding for the ATLAS experiment.

```{note}
While `DD4hep` offers a descriptive language with a dedicated extension mechanism
that can be used by Acts to interpret the underlying geometry hierarchy and and structure,
there is no such guarantee when having the already as built `TGeo` geometry in hand.
Therefor a dedicated Acts configuration file based on `json` can be provided that allows
to specify parsing restrictions for sub detectors. 
```


### Layer building

`Surface` object that are to be grouped on a layer should be put into a
`SurfaceArray` and provided to the layer. Certain helper tools exist to ease the
translation and create appropriate binning structure: The`SurfaceArrayCreator`
can create cylindrical, disc-like & planar layers, where the dimensions of the
layer are determined by parsing the provided surfaces. Additionally, an envelope
covering the surfaces can be chosen.

```{note}
There exist standard layer builders that are desinged to build cylindrical, disk like 
and planar layers and perform the ordering of the surfaces onto those layers. These
builders are called from the top level translation entry points from either `TGeo` 
or `DD4hep`.
```


### Volume building, packing, and gluing

The philosophy of the `TrackingGeometry` is a fully connective geometry setup,
i.e. `TrackingVolume` objects are either pure containers for other contained
`TrackingVolume` instances (where the contained volumes fully fill the space of
the container volume), or are fully attached via the boundary surface mechanism.
The boundary surfaces then act as portals from one `TrackingVolume` into the
next one along the trajectory.

The process to create a fully connected tracking geometry is called glueing.
Wherever possible, common boundary surfaces are *shared*, where this is not
possible, they are *attached*.

![GlueingBC](/figures/geometry/GlueingBC.png)
![GlueingABC](/figures/geometry/GlueingABC.png)
![NavigationABC](/figures/geometry/NavigationABC.png)

For cylindrical detector setups, a dedicated `CylinderVolumeBuilder` is
provided, which performs a variety of volume building, packing and gluing.

```{note}
For most cylindrical detectors, there exist automated glueing and geometry building 
modules that take care of the glueing process.
```

## TrackingGeometry building using a KDTree and a Proto Description

For cylindrical detectors there exist a generic tracking geometry building module,
based on KDTree and a proto description.

This building procedure uses a `ProtoDetector` description which provides a 
high level description of layers and container volumes, together with some 
binning and ordering information.
This proto description is then used to assign surfaces that are provided to the 
`KDTreeTrackingGeometryBuilder` using an internal query to the KD-tree structure.
