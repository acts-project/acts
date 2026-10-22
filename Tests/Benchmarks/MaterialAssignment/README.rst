Gen1 ODD material assignment benchmark
=====================================

This benchmark measures assignment candidate finding, including exact surface
intersections, allocations and path-length sorting. It calls the current Core
IntersectionMaterialAssigner as the baseline. It does not measure material
interaction assignment, accumulation, ROOT I/O or end-to-end mapping throughput.
Production material mapping is unchanged.

Reproduce
---------

Configure a Release build with ACTS_BUILD_BENCHMARKS=ON, then run from the
repository root::

    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DACTS_BUILD_BENCHMARKS=ON
    cmake --build build --target ActsBenchmarkMaterialAssignment
    build/bin/ActsBenchmarkMaterialAssignment \
      Tests/Benchmarks/MaterialAssignment/odd-gen1-material-surfaces.txt 100000 7

Use the usual dependency setup for your machine. The executable only needs Core;
DD4hep and ROOT are needed to regenerate the fixture, not to run the benchmark.
The second and third arguments are rays per sample and timed repetitions.
Before timing, each alternative must return the same ordered surface pointers,
intersection positions (within 1e-9 mm), and directions as the baseline.

The driver uses seed 20261015 and uniform phi and eta in [-4, 4]. Samples are:

* ip: origin (0, 0, 0).
* beamspot: uniform x,y in [-0.05, 0.05] mm and z in [-50, 50] mm.
* displaced: uniform x,y in [-1200, 1200] mm and z in [-4000, 4000] mm.

These are synthetic ray distributions on real ODD material surfaces, not a
recorded Geant4 material sample. Timings exclude ray generation, correctness
checks and index construction. Each method is warmed up once, execution order
rotates between repetitions, and the report gives median/min/max microseconds
per ray. The checksum consumes the results. Exact-test counts are measured
separately from timing. Construction times are single measurements.

Fixture provenance
------------------

odd-gen1-material-surfaces.txt contains the 107 surfaces returned by gen1 ODD
trackingGeometry().extractMaterialSurfaces(): 15 cylinders and 92 discs. It
preserves input order, geometry IDs, bounds and complete placements. It stores
no material payload because candidate finding does not use that payload.

Exported on 2026-10-15 using ODD commit
 d70556f33b1c36abdf101a07ad8cf57eb56fbdf1
and an installed ACTS Python build from the checkout at
 01f640f35d3cf6199152ea2ee45fd0f47178d925.
The baseline and prototypes were freshly compiled from this worktree, based on
 d66c0ec59523e5cc6bebf71de35e938ac3cc8570.
The installed Python examples initializer referenced an unavailable Helix binding;
for this export only, that unused propagator entry was removed in memory. No
geometry code or on-disk installed files were modified.

To regenerate with a consistent ACTS Python/DD4hep/ROOT installation::

    python Tests/Benchmarks/MaterialAssignment/export_odd.py surfaces.txt

The exporter explicitly selects gen3=False and uses the default ODD material
decorator. It rejects unsupported surface types instead of silently dropping them.

Prototype
---------

IndexedMaterialAssigner implements IAssignmentFinder in the benchmark directory.
It supports surface assignments only and offers three modes:

* Octree: existing makeOctree with depth 6 over conservative Cartesian AABBs.
  Traversal uses skip pointers and an inclusive slab predicate with cached
  reciprocals and explicit zero-direction handling, avoiding the existing
  ray-box predicate's boundary-parallel and zero-thickness issues.
* RZ: binary hierarchy of annular cylindrical bounds, split at the median
  centre along the larger r/z extent, traversed with skip indices.
* LinearRZ: the same cylindrical predicate applied to each surface without a tree.

For a ray p+t*d, first clip t>=0 to the bounding z interval. On that interval,
r(t)^2 is a convex quadratic. Its minimum is at the clamped closest transverse
approach; its maximum is at an endpoint, or infinity for an unbounded ray with
nonzero transverse direction. Reject when that radial range misses [rmin,rmax].
No quadratic roots or per-bound square roots are needed. For origins exactly
on the symmetry axis, r(t)=t*|d_xy| and the predicate becomes a 2D slab test.
Checking only radial values at both endpoints is unsafe for displaced rays
that first travel inward and then outward.

The filter retains the existing exact Surface::intersect call, closestForward
selection, strictly positive path-length check and final ordering. Candidate
indices are sorted back into input order to preserve the baseline's tie
behaviour. Bounds have twice s_onSurfaceTolerance padding. Only exactly coaxial
cylinders and radial discs are indexed; unsupported or tilted/translated-axis
surfaces are always tested. Phi sectors are conservatively enclosed in full
shells. Bounds belong to the construction geometry context and must be rebuilt
after alignment changes. This is a prototype, not a new production default.

The octree has a particular weakness for this fixture: its strict x/y/z octant
comparisons send centres on splitting planes into octant zero. All 107 AABBs
have x=y=0, so this build does not spatially subdivide them. Cartesian bounds
also include the large radial holes inside cylindrical/disc surfaces.

Results
-------

results-macos-arm64.txt is the captured 100000-ray, seven-repetition run on
macOS arm64, AppleClang 21.0.0.21000334, -O3 -DNDEBUG, Eigen 5.0.1.
Core was built locally; existing generated particle/propagator headers were
reused through ACTS_CODEGEN_PREBUILT_DIR (their generator cache keys matched).

Median microseconds per ray::

    Sample       Baseline   Octree    RZ tree   RZ linear
    ip           1.148      1.154     0.831     0.816
    beamspot     1.155      1.161     0.888     0.872
    displaced    0.957      0.702     0.470     0.531

RZ linear is about 1.41x faster for IP rays and 1.32x for beamspot rays;
RZ tree is about 2.04x faster for displaced rays. Both RZ modes reduce exact
intersections from 107 to about 20 per beamline ray and 7.3 per displaced ray.
The Cartesian octree is approximately tied with the baseline near the beamline.
These timings support using cylindrical bounds; with only 107 surfaces, a tree
is not consistently better than the simpler linear cylindrical filter.

Validation: all three alternatives matched the baseline on all 300000 random
rays and 2230 axial/transverse/tangent/near-parallel/boundary rays. The driver also
checks unsupported/non-coaxial fallback, empty input, and mixed indexed/fallback
input. No end-to-end mapping speedup is claimed. A production follow-up should
benchmark recorded material tracks, define alignment-cache invalidation and
volume handling, and add focused unit tests before exposing this assigner in Core.

Angular index and cylinder fast path experiment
----------------------------------------------

The driver now also measures an angular index for the fixed origin (0,0,0),
a brute-force control using the prototype's allocation/dispatch machinery,
that control with a cylinder fast path, and the angular index with the fast path.

Angular mode projects every conservative r-z bounding rectangle to an interval
in cot(theta), using the extrema of z/r at its four corners. Sorted interval
endpoints partition the angle domain; each open cell stores its complete
candidate list. A query uses binary search and selects that list, including
both neighbouring lists at an exact endpoint. This covers both cylinders and
discs, retaining input order and falling back to unconditional candidates for
radial bounds touching zero. Construction/storage can be quadratic in the
surface count because candidate lists overlap. This is reasonable for the 107
ODD surfaces but would need revisiting for a much larger detector.

The angular cache applies only to its fixed origin. All other origins, and
axial rays, use the linear r-z filter. Thus the beamspot/displaced angular rows
measure fallback performance, not angular acceleration at those origins.
Changing the fixed on-axis origin would require rebuilding z-relative intervals.

OnAxisCylinderIntersection is a separate compiled function, avoiding an unfair
inline-only comparison with the existing Core shared-library implementation.
For full cylinders with an exactly on-axis origin it computes both roots
+/-R/|direction cross axis|, then preserves the original on-surface, reachable
and unreachable statuses and z-boundary checks. Partial-phi cylinders, other
boundary tolerances, off-axis starts, axial directions and non-finite path
lengths use the existing intersection. Cylinder types are cached at index
construction, avoiding dynamic_cast per query. Neither production Core nor the
production mapper is modified.

The direct cylinder benchmark compares both intersections, their positions,
path lengths and statuses, with random on/off-axis starts, finite and infinite
bounds checks, another boundary tolerance, non-unit/axial/near-axial directions,
large on-surface tolerances, translated cylinders and partial-phi fallback.
Numerical comparisons allow roundoff; paths and positions need not be bitwise
identical to the quadratic solver. Timing uses separate on-axis and off-axis
samples with both statuses and the forward path consumed by a checksum.

The assignment driver additionally tests each angular interval endpoint and
its immediate floating-point neighbours (3514 edge rays in total), alongside
all 300000 random rays. Ordered assignment checks and fallback tests must pass
before the driver completes successfully.

Angular experiment results
--------------------------

results-angular-macos-arm64.txt captures the seven-repetition run with the same
fixture, compiler settings and sample definitions. Median IP assignment costs::

    Method                     microseconds/ray
    Core baseline              1.168
    Linear RZ                  0.837
    Angular index              0.621
    Angular + cylinder path    0.598
    Brute prototype control    1.236
    Brute + cylinder path      1.180

The angular index alone is 1.88x faster than Core and 1.35x faster than linear
RZ. The combined prototype is 1.95x faster than Core. The cylinder fast path
adds about 4% to angular performance; only 15/107 surfaces are cylinders. The
brute prototype control isolates allocation and dispatch overhead: its cylinder
variant improves by about 4.5%, but that does not beat the existing Core baseline.

The separate cylinder API experiment reports::

    Origins       Core ns/call   Prototype ns/call
    On axis       10.206         7.977
    Off axis      10.202         13.384

That is a 1.28x on-axis improvement, with a 31% off-axis regression from the
prototype's eligibility checks and fallback call. This wrapper is not a patch
to CylinderSurface::intersect; an integrated implementation could avoid some
fallback overhead and must be measured independently before a production change.
The direct test uses a fixed cylinder and None boundary tolerance for timing;
other geometry/tolerance cases are correctness checks, not performance claims.

For beamspot/displaced samples the angular index falls back to linear RZ; the
cylinder optimization offers no reliable advantage there. Prefer the angular
index for a known fixed on-axis origin, and retain a general r-z route otherwise.
Production adoption also needs a decision about cache lifetime, context changes,
input-origin constraints and the quadratic storage of precomputed candidate lists.

Unmeasured follow-up
--------------------

The angular index still returns about 20 candidates per IP ray. Exact surface
intersections, result allocations, restoration of input order and path-length
sorting remain. Their individual contribution has not been profiled, so the
current measurements do not establish which remaining cost dominates.

For fixed-origin, full coaxial cylinders and discs, an exact angular index could
identify hits and calculate positions directly, bypassing general surface
intersection calls. This has not been implemented or benchmarked here. A next
experiment should profile the remaining work, then compare direct hit generation
and allocation/sorting changes while preserving boundary tolerances, ordering
and empty-hit correction inputs.
