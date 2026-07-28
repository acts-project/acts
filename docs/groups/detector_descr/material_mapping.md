@defgroup material_mapping Material mapping
@brief Projection procedure to derive mapped material properties.
@ingroup material

Material mapping projects the material of a detailed simulation geometry (as
recorded by Geant4, see @ref material_mapping_workflow) onto the surfaces of a
tracking geometry, so that track reconstruction can account for material effects
without the detailed geometry being present.

The procedure is split into two independent steps, each behind its own
interface, and tied together by @ref Acts::MaterialMapper.

1. **Assignment** — @ref Acts::IAssignmentFinder decides, for a given recorded
   material track, which surfaces (and volumes) the individual material
   interactions along that track should be attributed to.
2. **Accumulation** — @ref Acts::ISurfaceMaterialAccumulator collects the
   assigned interactions per surface, bins them, and averages them into the
   final @ref Acts::ISurfaceMaterial. The default implementation is
   @ref Acts::BinnedSurfaceMaterialAccumulator.

Because the two steps only communicate through
@ref Acts::IAssignmentFinder::SurfaceAssignment, the way candidate surfaces are
found is a free choice, which is what the two available assigners differ in.

**Intersection-based assignment (default).**
@ref Acts::IntersectionMaterialAssigner is preconditioned with a flat list of
candidate surfaces and volumes, and for each material track it simply
intersects the track's ray with every candidate. It does not use the
@ref Acts::Navigator or a propagator at all, and correspondingly does not
depend on the geometry being navigable: the only thing it needs is a list of
material-carrying surfaces. This makes it agnostic to how the tracking geometry
was built — the same code path serves both the Gen1 and the Gen3 (blueprint)
construction, since both produce an @ref Acts::TrackingGeometry from which the
surface list can be extracted by a plain surface visitor. It also removes
navigation as a source of error: material that a navigator would have missed is
still assigned. The cost is that every candidate surface is tested for every
track, so the runtime scales with the number of candidates.

**Propagation-based assignment.**
@ref Acts::PropagatorMaterialAssigner runs a straight-line propagation through
the geometry and collects the surfaces the navigator actually encounters. This
is cheaper for large geometries, but it inherits whatever the navigation
configuration does or does not find, and it requires a navigable geometry.

The intersection-based assigner is the default in the Examples workflow. The
propagation-based one remains available and is useful as a cross-check.

**Deprecated interface.** @ref Acts::SurfaceMaterialMapper and
@ref Acts::VolumeMaterialMapper are the older, monolithic mappers that combined
assignment and accumulation and were hard-wired to a propagator and navigator.
They are deprecated in favour of @ref Acts::MaterialMapper composed with the
assigner and accumulator described above.

Related components:

- @ref Acts::MaterialInteractionAssignment implements the assignment logic that
  turns candidates into final assignments, including vetoes and re-assignments.
- @ref Acts::MaterialValidator re-records material from an already mapped
  geometry so that the mapped result can be compared against the Geant4 input.
- @ref material_mapping_workflow documents the end-to-end Examples chain
  (recording, mapping, validation) and the scripts that drive it.
