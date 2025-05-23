# Geometry alignment and misalignment

:::{note}
This section does not describe the process of *aligning* the detector, by the use of a track based alignment algorithms. It only describes how the contextual alignment is handled within the geometry description itself.
:::

Real life experiments often suffer from imperfections that have to be taken into account in order to achieve highest quality track reconstruction output. One of those imperfections is the fact that the actual detector elements might not be at the exactt location as they are thought to be, often as a result of installation uncertainty or the movement of detector components in time, e.g. through thermal expansion and contraction effects. The actual positioning of the detection elements is thus contextual data, i.e. it depends on the current (sometimes time-dependent) context.

## The {class}`Acts::GeometryContext` mechanism

In ACTS, this is handled with the {class}`Acts::GeometryContext` concept, which is a requirement for each geometric operation that relates detector elements (mostly expressed through a {class}`Acts::Surface` object).
As an example, the local to global transformation that relates a local measurement on a planar surface to its global 3D coordinates looks as follows:

:::{doxygenfunction} Acts::PlaneSurface::localToGlobal
:outline:
:::

The {class}`Acts::GeometryContext` object guarantees hereby that the geometric operation is executed within the contextual environment under which it is requested. Detector alignment, however, is a very specific aspect, and often bound not only to the primary geometry model (e.g. `GeoModel`, `DD4hep`, `ROOT::TGeo`, etc.), but also to the experiments conditions data structure and handling. ACTS thus offers a mechanism that allows to connect whatever infrastructure client code may have in place to deal with such geometric updates.

{class}`Acts::Surface` objects that represent detection elements are required to be connected to an object of a class extening the {class}`DetecctorElementBase` interface, and, if a surface is asked for its {class}`Acts::Transform3` objects that positions the surface within the global reference frame it simply relays this request to the corresponding detector element. Only at that stage the {class}`Acts::GeometryContext` which is based on `std::any` is unpacked and whatever information used to retrieve the correct contextual {class}`Acts::Transform3` object.

Client code can thus steer the behavior of this mechanism by implementing something along those lines (in this case the alignment context acts as a payload carrying the updated {class}`Acts::Transform3`):

```c++

class MyDetectorElement : public Acts::DetectorElementBase {

    public :
        /// This is a specific alignment context struct
        /// that carries the updated transform for this context
        struct AlignmentContext {
            /// Updated Transform set
            std::map<std::string, Acts::Transform3D> alignedTransforms;
        };

        /// Override: unpack the alignment context
        /// @param gctx the geometry context, will be unpacked if possbile
        const Acts::Transform3D& transform(const Acts::GeometryContext& gctx) const override;

    private:
        //!< The nominal transform
        Acts::Transform3 m_nominalTransform;
        //!< A unique name
        std::string m_name;

};
```

When a possible implementation of this function could look like:

```c++

const Acts::Transform3D& MyDetectorElement::transform(const Acts::GeometryContext& gctx) const {
    // Check if the geometry context has a value
    if (gctx.hasValue()){
        const AlignmentContext* myContext = gctx.maybeGet<AlignmentContext>();
        if (myContext != nullptr){
            auto fTransform = myContext->alignedTransforms.find(m_name);
            if (fTransform != myContext->alignedTransforms.end()){
                return fTransform->second;
            }
        }
    }
    return m_nominalTransform;

}
```

:::{note}
The lookup with `std::string` might possibly not be the most performant way to pull out the module's transform from the payload.
:::

In a very similar manner, the `AlignmentContext` could also connect to an external source for fetching the correct object.


## A possible implementation using the {class}`Acts::ITransformStore` and {class}`Acts::Delegate` classes

A possible implementation that can be used (and is indeed implemented in the detector elements in the `Plugin` directory) is given through a simple delegation mechanism:
the detector element extensions based on `DD4hep` (via `TGeo`), `GeoModel` and `Geant4` as defined in the corresponding plugins all implement this mechanism, in the case for `Geant4`, this looks as follows:

```c++
const Transform3& Geant4DetectorElement::transform(
    const GeometryContext& gctx) const {
  // This uses the contextual transform mechanism based on the proposed
  // AlignmentDelegate infastructure
  const Acts::Transform3* aTransform = Acts::contextualTransform(gctx, *this);
  if (aTransform != nullptr) {
    return *aTransform;
  }
  return m_toGlobal;
}
```

The unpacking of the {class}`Acts::GeometryContext` object is hereby happing in the `contextualTransform` helper function.


## Examples and showcases

The unit test `ActsUnitTestAlignmentContext` located in `Tests/UnitTests/Core/Geometry/AlignmentContextTests.cpp` showcases and tests both implementation, the extension of a detector element with a dedicated context class as shown above, and the usage of the generic `AlignmentDelegate` mechanism. They can, in fact, exist side by side without interference, as demonstrated.
