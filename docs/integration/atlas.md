# Using Acts from Athena

First steps have been taken towards use *Acts* inside *Athena*. The current goal
is to demonstrate that the Acts geometry design can accomodate current and
future ATLAS tracking geometries without requiring dedicated modifications in
the library.  This means building an `Acts::TrackingGeometry`, which contains a
hierarchy of `TrackingVolume`s. These, in turn, encompass layers of detector
elements.  Having done this, `Acts` can be used to propagate particles through
the geometry. This allows the testing of the navigation that is required to
follow a particle trajectory through the detector.

Since Acts is designed to be configurable, the construction of the layer
structure is delegated to a `LayerBuilder` utility object. A `LayerBuilder`
object, that is aware of Athena and the way the detector geometry is acessible
within it, can provide Acts with all the information that is needed to construct
the tracking geometry.

```cpp
Acts::SurfaceArrayCreator::Config sacCfg;
sacCfg.surfaceMatcher = matcher; // <- allows injection of Identifier info

auto surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>(
sacCfg,
Acts::getDefaultLogger("SurfaceArrayCreator", Acts::Logging::VERBOSE));

...

GeoModelLayerBuilder::Config cfg;
// layer builder gets the detector manager
cfg.mng = static_cast<const InDetDD::SiDetectorManager*>(manager);
gmLayerBuilder = std::make_shared<const GMLB>(cfg,
Acts::getDefaultLogger("GMLayBldr", Acts::Logging::VERBOSE));

Acts::CylinderVolumeBuilder::Config cvbConfig;
cvbConfig.layerEnvelopeR = {0, 0};
cvbConfig.layerEnvelopeZ       = 2;
cvbConfig.trackingVolumeHelper = cvh;
// ...
cvbConfig.layerBuilder         = gmLayerBuilder;

Acts::TrackingGeometryBuilder::Config tgbConfig;
tgbConfig.trackingVolumeHelper   = cylinderVolumeHelper;
tgbConfig.trackingVolumeBuilders = volumeBuilders;
auto trackingGeometryBuilder
  = std::make_shared<const Acts::TrackingGeometryBuilder>(tgbConfig);

std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry
  = trackingGeometryBuilder->trackingGeometry();
```

Tracking surfaces in Acts are coupled to a detector element object, which needs
to fulfill a certain interface. A dedicated detector element class, which can
wrap an ATLAS detector element object, was written. At this point, the class
tries to abstract away some differences between detector elements from the
different ATLAS subdetectors, which are irrelevant to Acts and allow the
`LayerBuilder` to be less verbose. The detector identification helpers are an
example for this. The detector element needs to provide a transform which
describes the transition into the local reference frame of it's active surface.
In this case, the transform is stored as a pointer and returned as a const
reference when requested. This should in principle allow for changes of the
transforms transparently to Acts. This would be helpful as an initial approach
to support handling of alignment, since no dedicated infrastructure for
alignment is in place within Acts.

## Current ATLAS Inner Detector

The detector managers allow access to all the detector elements. The layer
builder loops over all elements and wraps them into a Acts detector element.
Subsequently, the identifiers attached to the detector elements allow figuring
out which detector part and layer the elements belong to. The elements are thus
collected into buckets for each layer in each component. The elements are then
handed over to the *LayerCreator* utility, which is part of the Acts core, and
is versatile enough to handle ATLAS geometries. Identifier information can be
injected into the *LayerCreator*, which allows determining layer size and
binning automatically.

In the case of the TRT, the layer builder currently uses the `TRT_Numerology` to
figure out the number of layers without looping over elements. The detector
manager provides a method which allows access to detector elements based on the
indices returned by `TRT_Numerology`. As of now, every straw is translated into
a detector element in Acts. The endcaps are built as 160 separate disc layers
containing one layer of straws. In the barrel, there is one layer per barrel
ring. The binning is not optimal at this point, an arbitrary $\phi$ binning
should work reasonably well.

## ITk

With some modifications, the layer builder can build the ITk geometry in
`20.20`. One adjustment is the fact that the meaning of the identifier parts
*eta_module* and *layer_disk* changed for the Pixel endcaps. In the ID,
*layer_disk* enumerates the endcap themselves, while *eta_module* enumerates the
rings in $r$. In `20.20`, *eta_module* enumerates the $z$-position of rings,
while *layer_disk* enumerates the rings in $r$. A problem arises here, because
*eta_module* now only enumerates rings present in a given *eta_module*. This
way, the identifier cannot be used which elements belong to the same layer (in
the Acts sense).

## Design

Currently, the tracking geometry is built from inside one algorithm. The
particle propagation loop also runs within the algorithm. This will need to be
reworked, so that access to the Acts tracking geometry is available in a
reusable way.
