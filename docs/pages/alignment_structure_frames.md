@page alignment_structure_frames Frames for alignable structures

# Frames for alignable structures

`ActsAlignment::AlignableStructure` groups surfaces and defines which rigid-body
degrees of freedom may float. To use the composite Jacobians from
@ref alignment_composite_jacobians, it also needs an explicit frame: an origin
for the rotation pivot and an orientation for its local axes.

## Frame conventions

`localToGlobalTransform(gctx)` returns the structure-to-global rigid transform.
Its translation is the pivot in global coordinates; its rotation maps the
structure's local axes to global axes. Structure corrections use local
translations and local rotations. Surface alignment parameters still use global
translations and local rotations, so chaining to a surface derivative requires
`Acts::detail::localFrameToAlignmentParametersJacobian` as well as
`Acts::detail::compositeToComponentJacobian`.

Every structure frame is global, including child frames. `addChild()` declares
membership; it does not compose frames or move surfaces. The user chooses the
pivot explicitly. It is not computed from a centroid or inferred from a
geometry identifier.

## Fixed and contextual frames

A fixed frame is copied into the structure at construction:

```cpp
ActsAlignment::AlignableStructure structure(id, structureToGlobal);
```

It returns the same frame for every context. This is useful for a fixed geometry
or an explicitly chosen linearisation point. It does not follow later changes
to the supplied transform or to member surfaces.

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
