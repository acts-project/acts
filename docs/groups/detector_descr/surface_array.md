@defgroup surface_array Surface array
@ingroup geometry
@brief Binned lookup from a position on a layer to the sensitive surfaces near it

@ref Acts::SurfaceArray answers one question: given where a track crosses a
layer, which sensitive surfaces could it hit. It holds the surfaces of one
layer in a two-dimensional grid drawn on a representative surface and returns,
for a position and a direction, the surfaces registered in the bins around the
crossing point.

## Overview

A layer holds its sensitive surfaces in a `SurfaceArray`, and navigation asks
the array instead of testing every module:

- @ref Acts::SurfaceArrayCreator builds it, through
  @ref Acts::SurfaceArrayCreator::surfaceArrayOnCylinder,
  @ref Acts::SurfaceArrayCreator::surfaceArrayOnDisc or
  @ref Acts::SurfaceArrayCreator::surfaceArrayOnPlane, usually called from
  @ref Acts::LayerCreator.
- @ref Acts::Layer owns it and queries it from
  @ref Acts::Layer::compatibleSurfaces. This is the Gen1 navigation path.
- @ref Acts::SurfaceArrayNavigationPolicy builds one for a Gen3 layer volume
  and feeds the result into the navigation stream.

Both consumers call @ref Acts::SurfaceArray::neighbors with the geometry
context, the position and the direction, and receive a `std::span` of surface
pointers. The array does not intersect the surfaces it returns. It narrows the
candidates; the caller intersects them.

The single-surface constructor builds an array without a grid. Every lookup
returns that one surface.

## The grid

The grid is drawn on a representative @ref Acts::RegularSurface that the
creator derives from the extent of the surfaces:

| Layer    | Representative surface     | Axis 0        | Axis 1        |
|----------|----------------------------|---------------|---------------|
| Cylinder | @ref Acts::CylinderSurface | phi, `Closed` | z, `Bound`    |
| Disc     | @ref Acts::DiscSurface     | r, `Bound`    | phi, `Closed` |
| Plane    | @ref Acts::PlaneSurface    | first in-plane direction, `Bound` | second in-plane direction, `Bound` |

The axes are @ref Acts::IAxis objects, equidistant or variable, and every
lookup implementation is templated on both axis types behind the
`ISurfaceGridLookup` interface. The boundary type (@ref Acts::AxisBoundaryType)
fixes what happens at the edge. A `Closed` axis wraps: the bin after the last
phi bin is the first, and a neighborhood of an edge bin continues on the other
side. A `Bound` axis clamps: a position past the edge lands in the edge bin and
a neighborhood stops there. Each axis also carries an underflow and an
overflow bin, so @ref Acts::SurfaceArray::size returns `(n0 + 2) * (n1 + 2)`.

Bin edges are expressed in the local frame of the representative surface. On a
cylinder and a disc the creator rotates that frame so that the phi seam at
+-pi falls between modules, never on a module center.

### Surface local versus grid local

The grid is binned in the quantities the axes name. The bound-local
coordinates of the representative surface are not always those: a
@ref Acts::CylinderSurface measures arc length `R * phi` along its first local
coordinate while the grid bins `phi`. `surfaceToGridLocal` and
`gridToSurfaceLocal` divide and multiply by the cylinder radius; on a disc and
a plane they are the identity. Both maps are linear, so they carry a
displacement the same way as a position.

## Filling

The array registers every surface in every bin its projection onto the
representative surface overlaps. `fillSurfaceFootprint` does this once per
surface at construction:

1. Confirm that the surface belongs to the layer: its reference position must
   project onto the representative surface within the layer tolerance.
2. Take the outline from @ref Acts::Surface::polyhedronRepresentation and walk
   each edge in 32 samples. A straight edge in space is not straight in grid
   coordinates, so the vertices alone are not enough.
3. Project each sample along the normal of the representative surface onto
   it, convert to grid coordinates and to a bin. Samples in the underflow or
   overflow bins are dropped.
4. Group the sampled bins into columns keyed by the closed axis. Per column
   keep the lowest and highest bin along the other axis, and register the
   surface in every bin of every column span.

The span fill is exact when the projected outline is convex per column, which
holds for the module shapes in use. Keying the columns by the axis that wraps
means no span is ever reconstructed across the phi seam: a module straddling
+-pi produces columns at both ends of the phi axis, each with its own span
along the bound axis.

![A trapezoidal disc module projected onto the (phi, r) grid: the sampled outline, the bins the span fill registers, and a second module split at the phi seam.](geometry/surface_array_footprint.svg){width=600px}

After all surfaces are filled the bin contents are sorted and deduplicated,
`checkGrid` verifies that every input surface landed in at least one bin, and
the neighbor cache is built. The filling grid is scratch; its contents
survive as the zero-distance packs of the cache.

## Lookup

@ref Acts::SurfaceArray::neighbors "neighbors(gctx, position, direction)" is
the navigation entry point:

1. Intersect the representative surface along `direction` from `position`
   (`findCrossing`). No intersection means an empty result.
2. Convert the crossing to grid coordinates and to a bin.
3. Compute the neighbor distance per axis (`crossingNeighborDistance`).
4. Return the cached pack for that bin and distance.

@ref Acts::SurfaceArray::at "at(gctx, position, direction)" skips step 3 and
returns the content of the crossing bin alone.

### Why a window

The lookup happens where the track crosses the representative surface. A
module is registered where it projects onto that surface. The two differ: a
module sits at a radial offset from the representative surface, and between
the surface and the module the track moves along the layer. The shallower the
crossing, the further it moves. The window around the crossing bin spans that
movement.

![Lookup geometry in the r-z plane: staggered modules at two radii, the representative cylinder between them, an inclined track, and the slide in z between crossing the representative surface and reaching the module.](geometry/surface_array_lookup_rz.svg){width=640px}

`crossingNeighborDistance` derives the slide from the crossing geometry:

- The array is built with a `tolerance`, the half thickness of the layer. A
  track with direction `d` meeting the representative surface with normal `n`
  travels `tolerance / |n . d|` inside the layer on each side of the surface.
- Dropping the component along `n` leaves the part of that path that runs
  along the layer: `slide = (tolerance / |n . d|) * (d - (n . d) n)`.
- The slide is a global displacement.
  @ref Acts::Surface::localCartesianToBoundLocalDerivative maps it into the
  bound-local frame of the representative surface with the curvilinear scaling
  included: `(d(R phi), dz)` on a cylinder, `(dr, dphi)` on a disc,
  `(dx, dy)` on a plane. `surfaceToGridLocal` turns that into a step in grid
  coordinates.
- The bins at `+slide` and `-slide` from the crossing point give the distance
  per axis, in bins, with the wrap of a closed axis taken into account. The
  larger of the two sides is kept.

There is no case analysis on the surface type. The derivative supplies the
local metric for any @ref Acts::RegularSurface, and the result is already in
the quantities the grid is binned in.

![The same crossing on the (phi, z) grid: the crossing bin, the window opened along z only, the registered footprint of the module the track hits, and where the track meets it.](geometry/surface_array_lookup_grid.svg){width=640px}

Below an incidence `|n . d|` of `1e-4` the track runs along the layer and the
window opens to its bound.

## The neighbor window

@ref Acts::SurfaceArray::NeighborWindow bounds the distance per axis:

```cpp
struct NeighborWindow {
  std::array<std::uint8_t, 2> min = {0, 0};
  std::array<std::uint8_t, 2> max = {2, 2};
};
```

The computed distance is clamped into `[min, max]` on each axis. `max` also
sizes the neighbor cache, so it is a cost as well as a bound. The defaults set
by @ref Acts::SurfaceArrayCreator and @ref Acts::LayerCreator are:

| Layer    | `max`    | Meaning                      |
|----------|----------|------------------------------|
| Cylinder | `{1, 2}` | one bin in phi, two in z     |
| Disc     | `{2, 1}` | two bins in r, one in phi    |
| Plane    | `{2, 2}` | two bins in either direction |

`min` is `{0, 0}` everywhere. The floor exists for lookups under a geometry
context the fill did not see. The fill runs once, in the construction context;
a lookup under an alignment context sees moved surfaces, and a floor of one bin
keeps them reachable. `{1, 1}` for both `min` and `max` reproduces the fixed
3x3 lookup the array used before the window was sized by the crossing angle.
When `min` equals `max` the angle is not evaluated.

Constructing a `NeighborWindow` from a single `std::uint8_t` is deprecated; it
sets `max` to that value on both axes and `min` to one, or to zero if the bound
is zero.
@ref Acts::SurfaceArray::maxNeighborDistance is deprecated in favor of
@ref Acts::SurfaceArray::neighborWindow.

## The neighbor cache

The set of surfaces reachable from a bin at a given distance is fixed once the
array is filled. `populateNeighborCache` precomputes it for every valid bin and
every `(distance0, distance1)` up to `max`. Each result is a sorted,
deduplicated pack of surface pointers. Identical packs are stored once, found
through a hash table that reads packs back through their offsets, and the index
array holds `(offset, count)` pairs of 32-bit integers.

A lookup is therefore one index computation,
`bin * stridePerBin + distance0 * stride1 + distance1`, and a span into the
pack storage. No gather over neighboring bins happens at lookup time.

The index array has `(n0 + 2) * (n1 + 2) * (max0 + 1) * (max1 + 1)` entries
and is cheap. The pack contents are what cost memory: each pack lists the
surfaces of up to `(2 max0 + 1) * (2 max1 + 1)` bins. On the generic detector
at the default bounds the cache is 13.6 MB. A lookup on a toy barrel costs
about 80 ns.

## Evolution

- February 2018: `SurfaceArray` replaced the `BinUtility`/`BinnedArray`-based
  layer lookup. It already cached a neighbor list per bin.
- April 2018: boost type erasure removed from the array and from the axis
  introspection, giving the `ISurfaceGridLookup` pimpl with axis-templated
  implementations that is still in place.
- `SurfaceArray` is Gen1 geometry but not legacy: Gen3 reaches it through
  @ref Acts::SurfaceArrayNavigationPolicy.
- April 2026, PR #5186: the lookup became inclination dependent. The neighbor
  distance was derived from `1 / |n . d|` instead of being a fixed one bin.
- June 2026, PR #5498: deprecated members removed.
- July 2026, PR #5531: `MultiAxis` extracted from `Grid`; the lookup binds its
  axes through it.
- September 2026, PRs #6017, #6018 and #6019: the phi grid rotated into the
  frame of the representative surface, the axes built in that frame, and the
  fill changed from two point tests to the projected footprint.
- PR #6038: the neighbor window sized by the slide along the layer, per axis,
  bounded by @ref Acts::SurfaceArray::NeighborWindow.
