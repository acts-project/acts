@page alignment_structure_frames Frames for alignable structures

# Frames for alignable structures

`ActsAlignment::AlignableStructure` groups surfaces and defines which rigid-body
degrees of freedom may float. To use the composite Jacobians from
@ref alignment_composite_jacobians, it also needs an explicit frame: an origin
for the rotation pivot and an orientation for its local axes.

## Frame conventions

`localToGlobalTransform(gctx)` returns the structure-to-global rigid transform.
Its translation is the pivot in global coordinates; its rotation maps the
structure's local axes to global axes. Structure corrections and surface
alignment parameters both use local translations and local rotations, each in
their own frame. `Acts::detail::compositeToComponentJacobian` maps structure
corrections to the member surface's local parameters, including the rotation
lever arm. `Acts::detail::localFrameToAlignmentParametersJacobian` is identity
under this convention.

Every frame returned by `localToGlobalTransform(gctx)` maps structure-local
coordinates to global coordinates, including for child structures. This return
convention does not require applications to store or define geometry in global
coordinates. An application can express a child relative to its parent and
compose the transforms in the child's callable:

```cpp
const Acts::Transform3 childToGlobal = parentToGlobal * childToParent;
```

Both transforms must describe the geometry for the same `GeometryContext`.
The relative transform can also be recovered as
`parentToGlobal.inverse() * childToGlobal`. This allows parent-relative geometry
descriptions and reuse of symmetric layouts while retaining a common global
representation for the derivative calculations. Alignment corrections still
use each structure's own local axes.

`addChild()` declares membership; it does not perform this composition or move
surfaces. Automatic composition through the hierarchy, coordinated updates of
member surfaces, and hierarchy constraints are not implemented by this frame
API. The user chooses the pivot explicitly; it is not computed from a centroid
or inferred from a geometry identifier.

## Context-independent and contextual frames

A context-independent frame is copied into the structure at construction:

```cpp
ActsAlignment::AlignableStructure structure(id, structureToGlobal);
```

It returns the same frame for every context. This fixes the frame returned by
`localToGlobalTransform()`, not the alignment degrees of freedom, which are
controlled by `alignmentMask()`. The overload is useful for an explicitly
chosen pivot and orientation in derivative calculations or tests, and as an
initial linearisation frame for the planned structure solver. A fixed global
reference frame is one possible use, but is not required.

The stored frame does not follow later changes to the supplied transform or to
member surfaces. A solver using it as an initial frame must maintain subsequent
frame updates separately. The planned working-state mechanism is described
below; the frame API does not yet implement iterative structure alignment.

For a context-dependent frame, supply a callable instead. For example, an
application can put its own structure transforms in a geometry-context payload:

```cpp
struct StructureFrames {
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3> transforms;
};

ActsAlignment::AlignableStructure structure(
    id, [id](const Acts::GeometryContext& gctx) -> Acts::Transform3 {
      return gctx.get<StructureFrames>().transforms.at(id);
    });
```

The provider must return a rigid transform consistent with the member surfaces
evaluated in the same context. It is called on every query and must be safe for
concurrent reads. The structure owns the callable, but captured references and
external conditions data remain the application's lifetime responsibility.
Returned transforms are values, so callers cannot modify the provider's state
through a returned reference.

An iterative solver must apply a correction consistently to the structure frame
and to all member surfaces before the next linearisation. The frame provider is
a read interface; it does not implement these updates or mutate nominal geometry.
Contextual conditions storage remains outside the Alignment library.

For the planned iterative solver, providers seed working frames once per
alignment call. The solver will update those working frames together with the
member transforms; it must not reset them by rereading a fixed initial provider
at each iteration. Final frames must be made available to seed a later alignment
call. This working-state and update mechanism is not implemented by the frame
API.

## Metadata-only structures and Python

The existing constructor `AlignableStructure(id)` remains available for hierarchy
declarations. Such a structure has `hasTransform() == false`, and querying its
frame throws `std::logic_error`. There is no implicit identity frame. An empty
callable is rejected with `std::invalid_argument` at construction.

Python exposes the same choices:

```python
from acts.examples.alignment import AlignableStructure

fixed = AlignableStructure(id, transform=structure_to_global)
contextual = AlignableStructure(id, transformProvider=frame_provider)
assert fixed.hasTransform
frame = contextual.localToGlobalTransform(gctx)
```

Providing a frame does not yet enable numerical structure alignment. The native
solver still needs shared parameter indexing, Jacobian accumulation, and rigid
updates of the members. The frame API supplies the geometry contract for that
integration; hierarchy validation and standalone-module alignment retain their
existing behavior.
