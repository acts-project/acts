// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/TrackingGeometryLayout.hpp"

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace ActsFatras::Synthetic {

namespace {

/// One real surface, reduced to where it sits and how far it reaches.
struct SurfaceEntry {
  /// Centre of the surface
  double r{};
  double z{};
  /// Envelope of its footprint
  double rMin{};
  double rMax{};
  double zMin{};
  double zMax{};
};

/// A surface's material as it varies along the coordinate the surface extends
/// in: r for a disc, |z| for a cylinder. Equal bands over `[min, max]`.
struct SampledProfile {
  double min{};
  double max{};
  std::vector<Acts::MaterialSlab> bands;
};

/// Mean of a set of slabs, i.e. what a track meets crossing one of them.
/// @param slabs what to average over, empty giving vacuum
/// @param count what to divide by, which is the number of samples taken and not
///        the number that came back with anything: an empty patch of a surface
///        is material the track does not meet
/// @return the mean slab
Acts::MaterialSlab meanSlab(const std::vector<Acts::MaterialSlab>& slabs,
                            std::size_t count) {
  if (slabs.empty() || count == 0) {
    return Acts::MaterialSlab::Vacuum(0.f);
  }
  Acts::MaterialSlab mean = Acts::MaterialSlab::combineLayers(slabs);
  mean.scaleThickness(1.f / static_cast<float>(count));
  return mean;
}

/// What material a surface carries, band by band along its extent.
///
/// Sampled rather than taken at a point: mapped material is binned along the
/// surface, and the local origin of a disc is its beam hole, where there is
/// none.
///
/// @param surface the surface to measure
/// @param gctx the geometry context
/// @param bands how many bands to sample the surface into
/// @return the profile, or nothing where the geometry gives it no material
std::optional<SampledProfile> materialOf(const Acts::Surface& surface,
                                         const Acts::GeometryContext& gctx,
                                         std::size_t bands) {
  const Acts::ISurfaceMaterial* material = surface.surfaceMaterial();
  if (material == nullptr) {
    return std::nullopt;
  }

  const bool disc = surface.type() == Acts::Surface::Disc;
  const Acts::Polyhedron polyhedron = surface.polyhedronRepresentation(gctx, 1);
  Acts::Extent extent;
  extent.extend(polyhedron.vertices.begin(), polyhedron.vertices.end(),
                {Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisZ},
                /*applyEnv=*/false);
  const Acts::AxisDirection wanted =
      disc ? Acts::AxisDirection::AxisR : Acts::AxisDirection::AxisZ;
  auto range = extent.range(wanted);

  // Which local coordinate the material is binned along, asked rather than
  // assumed: the ODD's endcap material is binned (phi, r), so putting the
  // radius where a disc's own convention has it lands in the phi bin.
  std::size_t component = disc ? 0 : 1;
  const std::vector<Acts::AxisDirection> axes = material->localAxisDirections();
  for (std::size_t a = 0; a < axes.size() && a < 2; ++a) {
    if (axes[a] == wanted) {
      component = a;
    }
  }

  // The other local coordinate is the azimuthal one: rphi on a cylinder, phi on
  // a disc. It is sampled and averaged over -- a mapped surface is binned in
  // both and varies by a factor of six around itself, so one slice is not its
  // average -- where the first is kept.
  const std::size_t across = 1 - component;
  const double azimuthal =
      disc ? std::numbers::pi
           : std::numbers::pi * extent.range(Acts::AxisDirection::AxisR).max();

  // A cylinder is folded onto |z|, which is what the synthetic layout is
  // symmetric in and what `SurfaceMaterial::at` is given.
  const double lower = disc ? range.min() : 0.;
  const double upper =
      disc ? range.max() : std::max(std::abs(range.min()), range.max());
  if (!(upper > lower)) {
    return std::nullopt;
  }

  constexpr std::size_t samples = 16;
  SampledProfile profile{lower, upper, {}};
  profile.bands.reserve(bands);
  bool any = false;
  for (std::size_t i = 0; i < bands; ++i) {
    const double along = lower + (static_cast<double>(i) + 0.5) /
                                     static_cast<double>(bands) *
                                     (upper - lower);
    // both ends of a cylinder, which fold onto the same |z|
    const std::array<double, 2> sides{along, -along};
    const std::size_t numSides = disc ? 1 : 2;
    std::vector<Acts::MaterialSlab> slabs;
    std::size_t taken = 0;
    for (std::size_t s = 0; s < numSides; ++s) {
      const double side = sides[s];
      if (side < range.min() || side > range.max()) {
        continue;
      }
      for (std::size_t j = 0; j < samples; ++j) {
        Acts::Vector2 local = Acts::Vector2::Zero();
        local[component] = side;
        local[across] = -azimuthal + 2. * azimuthal *
                                         (static_cast<double>(j) + 0.5) /
                                         static_cast<double>(samples);
        ++taken;
        const Acts::MaterialSlab& slab = material->materialSlab(local);
        if (!slab.isVacuum()) {
          slabs.push_back(slab);
        }
      }
    }
    any = any || !slabs.empty();
    profile.bands.push_back(meanSlab(slabs, taken));
  }
  return any ? std::optional{profile} : std::nullopt;
}

/// Put a profile onto a different set of bands, keeping what it integrates to.
/// @param profile the profile to resample
/// @param min the lower bound of the target bands
/// @param max the upper bound
/// @param bands how many equal bands to spread it over
/// @return the resampled bands, vacuum where the profile does not reach
std::vector<Acts::MaterialSlab> resample(const SampledProfile& profile,
                                         double min, double max,
                                         std::size_t bands) {
  const double width = (max - min) / static_cast<double>(bands);
  const double source =
      (profile.max - profile.min) / static_cast<double>(profile.bands.size());
  std::vector<Acts::MaterialSlab> out;
  out.reserve(bands);
  for (std::size_t i = 0; i < bands; ++i) {
    const double lo = min + static_cast<double>(i) * width;
    const double hi = lo + width;
    std::vector<Acts::MaterialSlab> parts;
    for (std::size_t k = 0; k < profile.bands.size(); ++k) {
      const double overlap =
          std::min(hi, profile.min + static_cast<double>(k + 1) * source) -
          std::max(lo, profile.min + static_cast<double>(k) * source);
      if (overlap <= 0. || profile.bands[k].isVacuum()) {
        continue;
      }
      Acts::MaterialSlab part = profile.bands[k];
      part.scaleThickness(static_cast<float>(overlap / width));
      parts.push_back(part);
    }
    out.push_back(parts.empty() ? Acts::MaterialSlab::Vacuum(0.f)
                                : Acts::MaterialSlab::combineLayers(parts));
  }
  return out;
}

/// A material surface together with where it sits, which is what says whether
/// two of them are the same crossing or two.
struct PositionedSlab {
  double r{};
  double z{};
  SampledProfile profile;
};

/// A group of surfaces that will become one or more synthetic surfaces.
struct SurfaceGroup {
  /// The r and z envelope of the group, which decides whether it is a cylinder
  /// or a disc and gives a cylinder its bounds
  Acts::Extent extent;
  bool sensitive{};
  /// The individual surfaces, which is what the ring structure of a disc is
  /// recovered from
  std::vector<SurfaceEntry> entries;
  /// The material of the layer, which a mapped geometry carries on its
  /// representing and approach surfaces. Those at different positions add up;
  /// coincident ones do not, and a group holds several of both.
  std::vector<PositionedSlab> layerMaterials;
  /// The material of the sensors, which a geometry built with it rather than
  /// mapped onto it carries there instead. A track crosses *one* of them, so
  /// these are averaged.
  std::vector<Acts::MaterialSlab> sensorMaterials;
};

/// The material of a layer: what its own surfaces add up to, plus what one of
/// its sensors carries, band by band along the layer.
/// @param group the group to sum over
/// @param bands how many bands the result carries
/// @return the profile, empty if the group has no material
SampledProfile materialOfGroup(const SurfaceGroup& group, std::size_t bands) {
  // A track crosses a *position* once. The two approach surfaces straddling a
  // layer's sensors are two crossings and add; the same boundary seen from the
  // volumes on either side is one, and adding it counts the material twice.
  // The ITk keeps six coincident cylinders at r = 124, which summed is most of
  // the pixel radiation length.
  constexpr double coincident = 1.;
  std::vector<PositionedSlab> positions;
  std::vector<std::vector<const SampledProfile*>> here;
  for (const PositionedSlab& entry : group.layerMaterials) {
    std::size_t at = positions.size();
    for (std::size_t i = 0; i < positions.size(); ++i) {
      if (std::abs(positions[i].r - entry.r) < coincident &&
          std::abs(positions[i].z - entry.z) < coincident) {
        at = i;
        break;
      }
    }
    if (at == positions.size()) {
      positions.push_back(entry);
      here.emplace_back();
    }
    here[at].push_back(&entry.profile);
  }

  // The surfaces of one layer do not all span the same range, an approach
  // surface reaching past the sensors it stands for, so they are summed over
  // the range that covers them all.
  SampledProfile out;
  out.min = std::numeric_limits<double>::max();
  out.max = std::numeric_limits<double>::lowest();
  for (const PositionedSlab& entry : group.layerMaterials) {
    out.min = std::min(out.min, entry.profile.min);
    out.max = std::max(out.max, entry.profile.max);
  }
  if (!(out.max > out.min)) {
    if (group.sensorMaterials.empty()) {
      return {};
    }
    // sensors carry no profile of their own, so they need a range to sit on
    out.min = 0.;
    out.max = 1.;
  }

  std::vector<std::vector<Acts::MaterialSlab>> stacked(bands);
  for (const std::vector<const SampledProfile*>& at : here) {
    // coincident surfaces are one crossing, so they are averaged before the
    // positions along the layer are added up
    std::vector<std::vector<Acts::MaterialSlab>> mean(bands);
    for (const SampledProfile* profile : at) {
      const std::vector<Acts::MaterialSlab> onto =
          resample(*profile, out.min, out.max, bands);
      for (std::size_t i = 0; i < bands; ++i) {
        if (!onto[i].isVacuum()) {
          mean[i].push_back(onto[i]);
        }
      }
    }
    for (std::size_t i = 0; i < bands; ++i) {
      const Acts::MaterialSlab slab = meanSlab(mean[i], at.size());
      if (!slab.isVacuum()) {
        stacked[i].push_back(slab);
      }
    }
  }
  if (!group.sensorMaterials.empty()) {
    // A track crosses one sensor of a layer, not all of them, and a sensor
    // spans too little of the layer to say where along it the material sits.
    const Acts::MaterialSlab sensor =
        meanSlab(group.sensorMaterials, group.sensorMaterials.size());
    for (std::size_t i = 0; i < bands; ++i) {
      stacked[i].push_back(sensor);
    }
  }

  out.bands.reserve(bands);
  bool any = false;
  for (std::size_t i = 0; i < bands; ++i) {
    out.bands.push_back(stacked[i].empty()
                            ? Acts::MaterialSlab::Vacuum(0.f)
                            : Acts::MaterialSlab::combineLayers(stacked[i]));
    any = any || !out.bands.back().isVacuum();
  }
  return any ? out : SampledProfile{};
}

/// Turn a sampled profile into what a surface carries: the composition and
/// thickness the bands share, and the two lengths of each of them.
///
/// Bands are merged so that the shipped tables stay readable, most of a real
/// surface being uniform. Both lengths have to agree for two to be one, so a
/// band holding as much material as its neighbour but of a different kind
/// survives.
///
/// @param profile the sampled profile
/// @param tolerance relative difference below which two bands are merged
/// @return the material
SurfaceMaterial compress(const SampledProfile& profile, double tolerance) {
  if (profile.bands.empty()) {
    return {};
  }
  const Acts::MaterialSlab mean = meanSlab(profile.bands, profile.bands.size());
  if (mean.isVacuum()) {
    return {};
  }

  const double width =
      (profile.max - profile.min) / static_cast<double>(profile.bands.size());
  // The first edge is where the material starts, so a surface whose material
  // begins away from the beam axis needs no band of nothing in front of it.
  std::vector<float> bounds{static_cast<float>(profile.min)};
  std::vector<float> radiationLengths;
  std::vector<float> nuclearLengths;

  /// What a run of sampled bands comes to, per unit of the shared thickness.
  /// @param from the first of the run
  /// @param to one past the last
  /// @return x/X0 and x/L0
  auto over = [&](std::size_t from, std::size_t to) {
    double x0 = 0.;
    double l0 = 0.;
    for (std::size_t i = from; i < to; ++i) {
      x0 += profile.bands[i].thicknessInX0();
      l0 += profile.bands[i].thicknessInL0();
    }
    const double n = static_cast<double>(to - from);
    return std::pair{x0 / n, l0 / n};
  };

  std::size_t start = 0;
  for (std::size_t k = 1; k <= profile.bands.size(); ++k) {
    const auto [x0, l0] = over(start, k);
    if (k < profile.bands.size()) {
      const auto [nextX0, nextL0] = over(k, k + 1);
      const auto apart = [tolerance](double a, double b) {
        const double hi = std::max(a, b);
        return hi > 0. && (hi - std::min(a, b)) > tolerance * hi;
      };
      if (!apart(x0, nextX0) && !apart(l0, nextL0)) {
        continue;
      }
    }
    bounds.push_back(
        static_cast<float>(profile.min + static_cast<double>(k) * width));
    // the band is quoted at the shared thickness, so its lengths are what give
    // that thickness the x/X0 and x/L0 the samples came to
    radiationLengths.push_back(
        x0 > 0. ? static_cast<float>(mean.thickness() / x0) : 0.f);
    nuclearLengths.push_back(l0 > 0. ? static_cast<float>(mean.thickness() / l0)
                                     : 0.f);
    start = k;
  }
  const Acts::Material& stuff = mean.material();
  const BandComposition composition{
      static_cast<float>(stuff.Ar()), static_cast<float>(stuff.Z()),
      static_cast<float>(stuff.molarDensity() * stuff.X0()),
      static_cast<float>(mean.thickness())};
  return SurfaceMaterial{composition, std::move(bounds), radiationLengths,
                         nuclearLengths};
}

/// Key a group by the volume and layer of the geometry identifier, which is the
/// grouping the tracking geometry already provides.
/// What makes two surfaces one layer: the subsystem they were assigned to, and
/// the volume and layer their geometry identifier puts them in.
struct GroupKey {
  std::string subsystem;
  Acts::GeometryIdentifier::Value volume{};
  Acts::GeometryIdentifier::Value layer{};

  /// @param other the key to compare against
  /// @return whether this one sorts before it
  bool operator<(const GroupKey& other) const {
    return std::tie(subsystem, volume, layer) <
           std::tie(other.subsystem, other.volume, other.layer);
  }
};

/// @param surface the surface to place
/// @param options the steering, for the subsystem it names
/// @return the group it belongs to
GroupKey groupKeyOf(const Acts::Surface& surface,
                    const TrackingGeometryLayoutOptions& options) {
  const Acts::GeometryIdentifier id = surface.geometryId();
  return {options.subsystemName ? options.subsystemName(surface)
                                : defaultSubsystemName(surface),
          id.volume(), id.layer()};
}

/// Where a surface sits in r.
///
/// A cylinder or a disc is centred on the beam axis and has no radius there, so
/// its envelope gives one. Anything else is a module and its centre says where
/// it is -- which is *inside* its own corner radii, a plane tangent to a
/// cylinder being nearer the beam at its middle than at its edges. Reading the
/// envelope instead puts the ODD's innermost barrel layer a millimetre out.
///
/// @param surface the surface
/// @param centre its centre
/// @param rMin the inner edge of its footprint
/// @param rMax the outer edge
/// @return the radius to place it at
double representativeR(const Acts::Surface& surface,
                       const Acts::Vector3& centre, double rMin, double rMax) {
  const bool aboutTheBeam = surface.type() == Acts::Surface::Cylinder ||
                            surface.type() == Acts::Surface::Disc;
  return aboutTheBeam ? 0.5 * (rMin + rMax) : Acts::VectorHelpers::perp(centre);
}

/// Accumulate a surface's footprint and centre into a group.
void extendBy(SurfaceGroup& group, const Acts::Surface& surface,
              const Acts::GeometryContext& gctx, std::size_t bands) {
  const std::vector<Acts::AxisDirection> axes{Acts::AxisDirection::AxisR,
                                              Acts::AxisDirection::AxisZ};
  // one segment per quarter is enough: we only want the r and z envelope, and
  // this keeps the vertex count down over a detector with 10^5 sensors
  const Acts::Polyhedron polyhedron = surface.polyhedronRepresentation(gctx, 1);
  Acts::Extent own;
  own.extend(polyhedron.vertices.begin(), polyhedron.vertices.end(), axes,
             /*applyEnv=*/false);
  group.extent.extend(own, axes, /*applyEnv=*/false);

  const Acts::Vector3 centre = surface.center(gctx);
  const double rMin = own.range(Acts::AxisDirection::AxisR).min();
  const double rMax = own.range(Acts::AxisDirection::AxisR).max();
  const double centreR = representativeR(surface, centre, rMin, rMax);
  group.entries.push_back(
      SurfaceEntry{centreR, centre.z(), rMin, rMax,
                   own.range(Acts::AxisDirection::AxisZ).min(),
                   own.range(Acts::AxisDirection::AxisZ).max()});
  if (const std::optional<SampledProfile> profile =
          materialOf(surface, gctx, bands);
      profile.has_value()) {
    if (surface.geometryId().sensitive() != 0) {
      group.sensorMaterials.push_back(
          meanSlab(profile->bands, profile->bands.size()));
    } else {
      group.layerMaterials.push_back(
          PositionedSlab{centreR, centre.z(), *profile});
    }
  }
}

/// Mean of a projection over the surfaces of a group. The synthetic surface
/// goes at the mean of the centres rather than the middle of the envelope,
/// which the tilt and corners of the real sensors pull off.
/// @param entries the surfaces to average over, which must not be empty
/// @param value what to average
/// @return the mean
template <typename Projection>
double meanOf(const std::vector<SurfaceEntry>& entries, Projection value) {
  double sum = 0.;
  for (const SurfaceEntry& entry : entries) {
    sum += value(entry);
  }
  return sum / static_cast<double>(entries.size());
}

/// What one group is turned into: a cylinder with the z range it spans, or a
/// disc with the rings it carries.
struct ReducedSurface {
  /// The subsystem it belongs to, which is what a layer of it is named by
  std::string subsystem;
  SurfaceShape shape{};
  /// signed z for a disc, radius for a cylinder
  double refCoord{};
  /// z range of a cylinder; unused for a disc
  double minBound{};
  double maxBound{};
  /// radial bands of a disc; unused for a cylinder
  std::vector<RingBounds> rings;
  /// The material of the layer, shared out over the surfaces reduced from it
  SurfaceMaterial material;
  bool sensitive{};
};

/// Group the surfaces of a disc into the radial bands they form: modules of one
/// ring overlap slightly in r and merge, the next ring out does not.
/// @param entries the surfaces of one z plane
/// @param tolerance the radial gap below which two bands are one
/// @return the bands, increasing in r
std::vector<RingBounds> ringsOf(std::vector<SurfaceEntry> entries,
                                double tolerance) {
  std::ranges::sort(entries, {}, &SurfaceEntry::rMin);

  std::vector<RingBounds> rings;
  for (const SurfaceEntry& entry : entries) {
    if (!rings.empty() &&
        entry.rMin - static_cast<double>(rings.back().rMax) <= tolerance) {
      rings.back().rMax =
          std::max(rings.back().rMax, static_cast<float>(entry.rMax));
    } else {
      rings.push_back(RingBounds{static_cast<float>(entry.rMin),
                                 static_cast<float>(entry.rMax)});
    }
  }
  return rings;
}

/// How much two sets of radial bands cover in common, as a fraction of the
/// narrower band. Zero says they are different rings of a disc, one says they
/// are the same ring twice.
/// @param a the bands of one z plane
/// @param b the bands of another
/// @return the largest overlap fraction over all pairs of bands
double overlapFraction(const std::vector<RingBounds>& a,
                       const std::vector<RingBounds>& b) {
  double worst = 0.;
  for (const RingBounds& ra : a) {
    for (const RingBounds& rb : b) {
      const double overlap =
          std::min(ra.rMax, rb.rMax) - std::max(ra.rMin, rb.rMin);
      if (overlap <= 0.) {
        continue;
      }
      const double narrower = std::min(ra.rMax - ra.rMin, rb.rMax - rb.rMin);
      if (narrower > 0.) {
        worst = std::max(worst, overlap / narrower);
      }
    }
  }
  return worst;
}

/// Reduce one group to the synthetic surfaces that stand for it: a cylinder is
/// always one, a disc one per z plane its modules fall into.
/// @param group the group to reduce
/// @param options steering
/// @return the surfaces
std::vector<ReducedSurface> reduce(
    const SurfaceGroup& group, const TrackingGeometryLayoutOptions& options) {
  const auto rRange = group.extent.range(Acts::AxisDirection::AxisR);
  const auto zRange = group.extent.range(Acts::AxisDirection::AxisZ);

  ReducedSurface reduced;
  reduced.sensitive = group.sensitive;
  reduced.material = compress(materialOfGroup(group, options.materialBands),
                              options.materialBandTolerance);

  // A barrel layer is thin in r and long in z, an endcap disc the other way
  // round. `Surface::type` is no help here: the sensitive surfaces of a real
  // tracker are planes in both cases. Both tests are needed: the ITk's
  // innermost endcap ring is twenty millimetres wide, thin enough in r to pass
  // for a barrel layer on its own, and only its z extent says otherwise.
  if (rRange.max() - rRange.min() < options.cylinderRTolerance &&
      rRange.max() - rRange.min() < zRange.max() - zRange.min()) {
    reduced.shape = SurfaceShape::Cylinder;
    reduced.refCoord =
        meanOf(group.entries, [](const SurfaceEntry& e) { return e.r; });
    reduced.minBound = zRange.min();
    reduced.maxBound = zRange.max();
    return {reduced};
  }

  reduced.shape = SurfaceShape::Disc;

  // split the layer into the z planes its modules sit on, each with the radial
  // bands those modules form
  std::vector<SurfaceEntry> entries = group.entries;
  std::ranges::sort(entries, {}, &SurfaceEntry::z);

  struct Plane {
    std::vector<SurfaceEntry> entries;
    std::vector<RingBounds> rings;
  };
  std::vector<Plane> planes;
  for (const SurfaceEntry& entry : entries) {
    if (!planes.empty() &&
        entry.z - planes.back().entries.back().z <= options.discZTolerance) {
      planes.back().entries.push_back(entry);
    } else {
      planes.push_back(Plane{{entry}, {}});
    }
  }
  for (Plane& plane : planes) {
    plane.rings = ringsOf(plane.entries, options.ringRTolerance);
  }

  // Merge back any two planes that cover the same radii: that is one ring
  // staggered along z, and splitting it would give a track two space points
  // where the detector gives it one. Every pair has to be tested, not just
  // neighbours in z, the generic detector's strip endcap interleaving its
  // rings.
  for (bool merged = true; merged && planes.size() > 1;) {
    merged = false;
    for (std::size_t p = 0; p + 1 < planes.size() && !merged; ++p) {
      for (std::size_t q = p + 1; q < planes.size(); ++q) {
        if (overlapFraction(planes[p].rings, planes[q].rings) <=
            options.maxRingOverlap) {
          continue;
        }
        planes[p].entries.insert(planes[p].entries.end(),
                                 planes[q].entries.begin(),
                                 planes[q].entries.end());
        planes.erase(planes.begin() + static_cast<std::ptrdiff_t>(q));
        planes[p].rings = ringsOf(planes[p].entries, options.ringRTolerance);
        merged = true;
        break;
      }
    }
  }

  std::vector<ReducedSurface> discs;
  discs.reserve(planes.size());
  for (const Plane& plane : planes) {
    // a disc whose surfaces straddle z = 0 averages to nothing, which the
    // caller filters out below
    reduced.refCoord =
        meanOf(plane.entries, [](const SurfaceEntry& e) { return e.z; });
    reduced.rings = subdivideRings(plane.rings, options.modulesPerSurface);
    // the radial extent, which is what bounds it when it has no rings
    reduced.minBound = rRange.min();
    reduced.maxBound = rRange.max();
    discs.push_back(reduced);
  }
  // the layer's material is shared out over the planes it was split into: a
  // track crosses the layer once whichever plane records it
  if (!discs.empty() && !discs.front().material.average().isVacuum()) {
    for (ReducedSurface& disc : discs) {
      disc.material.scaleThickness(1.f / static_cast<float>(discs.size()));
    }
  }
  return discs;
}

}  // namespace

DetectorDescription makeDescriptionFromTrackingGeometry(
    const Acts::TrackingGeometry& trackingGeometry,
    const Acts::GeometryContext& gctx,
    const TrackingGeometryLayoutOptions& options) {
  std::map<GroupKey, SurfaceGroup> groups;

  trackingGeometry.visitSurfaces(
      [&](const Acts::Surface* surface) {
        if (surface == nullptr) {
          return;
        }
        if (options.surfaceSelector && !options.surfaceSelector(*surface)) {
          return;
        }
        SurfaceGroup& group = groups[groupKeyOf(*surface, options)];
        group.sensitive = true;
        extendBy(group, *surface, gctx, options.materialBands);
      },
      /*restrictToSensitives=*/true);

  if (options.materialFromGeometry || options.includeMaterialSurfaces) {
    // Second pass over the surfaces that carry material, which in a mapped
    // geometry are the representing and approach surfaces of a layer and not
    // its sensors.
    trackingGeometry.visitSurfaces(
        [&](const Acts::Surface* surface) {
          if (surface == nullptr) {
            return;
          }
          if (surface->geometryId().sensitive() != 0) {
            // the first pass owns the sensors, whose material is averaged and
            // not summed: a track crosses one of them, not all
            return;
          }
          const std::optional<SampledProfile> profile =
              materialOf(*surface, gctx, options.materialBands);
          if (!profile.has_value()) {
            return;
          }
          // the selector applies here too, otherwise a layout restricted to one
          // subdetector still collects the material of the whole geometry
          if (options.surfaceSelector && !options.surfaceSelector(*surface)) {
            return;
          }
          const GroupKey key = groupKeyOf(*surface, options);
          const auto it = groups.find(key);
          if (it != groups.end()) {
            // Already covered by a sensitive group, which keeps its own
            // footprint: a layer surface spans more than the sensors it stands
            // for. Its material is taken, being the layer's.
            const Acts::Vector3 centre = surface->center(gctx);
            const Acts::Polyhedron shape =
                surface->polyhedronRepresentation(gctx, 1);
            Acts::Extent own;
            own.extend(shape.vertices.begin(), shape.vertices.end(),
                       {Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisZ},
                       /*applyEnv=*/false);
            it->second.layerMaterials.push_back(PositionedSlab{
                representativeR(*surface, centre,
                                own.range(Acts::AxisDirection::AxisR).min(),
                                own.range(Acts::AxisDirection::AxisR).max()),
                centre.z(), *profile});
            return;
          }
          if (!options.includeMaterialSurfaces) {
            // material away from any sensitive layer, which is a crossing of
            // its own and only wanted when asked for
            return;
          }
          SurfaceGroup& group = groups[key];
          group.sensitive = false;
          extendBy(group, *surface, gctx, options.materialBands);
        },
        /*restrictToSensitives=*/false);
  }

  // Reduce each group, then sort so that the builder's implicit numbering
  // (cylinders and discs in call order) comes out ordered by distance from the
  // interaction point.
  std::vector<ReducedSurface> cylinders;
  std::vector<ReducedSurface> discs;
  for (const auto& [key, group] : groups) {
    if (group.entries.empty()) {
      continue;
    }
    if (!group.extent.constrains(Acts::AxisDirection::AxisR) ||
        !group.extent.constrains(Acts::AxisDirection::AxisZ)) {
      // nothing was accumulated, e.g. a surface with no polyhedron
      continue;
    }
    for (ReducedSurface& reduced : reduce(group, options)) {
      // what the group was keyed by, which is what its layers are named by
      reduced.subsystem = key.subsystem;
      if (reduced.shape == SurfaceShape::Cylinder) {
        cylinders.push_back(std::move(reduced));
      } else {
        discs.push_back(std::move(reduced));
      }
    }
  }

  std::ranges::sort(cylinders, {}, &ReducedSurface::refCoord);
  // discs are numbered per side, counting outwards from the interaction point
  std::ranges::sort(discs,
                    [](const ReducedSurface& a, const ReducedSurface& b) {
                      return std::abs(a.refCoord) < std::abs(b.refCoord);
                    });

  DetectorDescription description;
  if (options.escapeRadius.has_value()) {
    description.escapeRadius = *options.escapeRadius;
  }
  if (options.escapeHalfZ.has_value()) {
    description.escapeHalfZ = *options.escapeHalfZ;
  }

  if (options.passiveBeamPipeRadius.has_value()) {
    // The beam pipe belongs to the detector rather than to any of its
    // subsystems, as it does in the geometry: it sits in a volume of its own,
    // which is why a subsystem selector excludes it and this option exists.
    description.passives.push_back(
        PassiveSurfaceDescription{.shape = SurfaceShape::Cylinder,
                                  .refCoord = *options.passiveBeamPipeRadius,
                                  .maxBound = description.escapeHalfZ});
  }

  /// The subsystem a reduced surface belongs to, made if it is new. Kept in the
  /// order the surfaces were found, so that a detector comes out ordered by
  /// distance from the interaction point rather than by name.
  const auto subsystemFor =
      [&description](const std::string& name) -> SubsystemDescription& {
    const auto existing = std::ranges::find(description.subsystems, name,
                                            &SubsystemDescription::name);
    if (existing != description.subsystems.end()) {
      return *existing;
    }
    description.subsystems.push_back(SubsystemDescription{.name = name});
    return description.subsystems.back();
  };

  /// Where a subsystem keeps the endcap of one side, made if it is new. One per
  /// side rather than one mirrored: a real detector is only nearly symmetric
  /// and this is a measurement of it, so what was found on either side is what
  /// is written down.
  const auto endcapFor = [](SubsystemDescription& subsystem,
                            const SurfaceSide side) -> EndcapDescription& {
    const EndcapPlacement placement = side == SurfaceSide::Positive
                                          ? EndcapPlacement::Positive
                                          : EndcapPlacement::Negative;
    const auto existing = std::ranges::find(subsystem.endcaps, placement,
                                            &EndcapDescription::placement);
    if (existing != subsystem.endcaps.end()) {
      return *existing;
    }
    subsystem.endcaps.push_back(EndcapDescription{.placement = placement});
    return subsystem.endcaps.back();
  };

  const auto materialOrNothing = [&options](const ReducedSurface& reduced) {
    return options.materialFromGeometry ? reduced.material : SurfaceMaterial{};
  };

  for (const ReducedSurface& cylinder : cylinders) {
    SubsystemDescription& subsystem = subsystemFor(cylinder.subsystem);
    const float radius = static_cast<float>(cylinder.refCoord);
    // the synthetic cylinder is symmetric in z, so it has to cover whichever
    // end of the real one reaches further
    const float halfLengthZ = static_cast<float>(
        std::max(std::abs(cylinder.minBound), std::abs(cylinder.maxBound)));
    if (!cylinder.sensitive) {
      // its own z extent, not the default of infinity: a service tube grazed
      // at cosh(eta) is worth ten times its thickness forward, and the real
      // one ends with the detector
      subsystem.passives.push_back(
          PassiveSurfaceDescription{.shape = SurfaceShape::Cylinder,
                                    .refCoord = radius,
                                    .maxBound = halfLengthZ,
                                    .material = materialOrNothing(cylinder)});
      continue;
    }
    if (subsystem.barrels.empty()) {
      subsystem.barrels.emplace_back();
    }
    subsystem.barrels.front().cylinders.push_back(
        CylinderDescription{.radius = radius,
                            .halfLengthZ = halfLengthZ,
                            .modules = options.modulesPerSurface,
                            .material = materialOrNothing(cylinder)});
  }

  for (const ReducedSurface& disc : discs) {
    // a disc centred on z = 0 belongs to no side, which a real endcap does not
    // produce
    if (disc.refCoord == 0.) {
      continue;
    }
    SubsystemDescription& subsystem = subsystemFor(disc.subsystem);
    const SurfaceSide side =
        disc.refCoord > 0. ? SurfaceSide::Positive : SurfaceSide::Negative;
    const float absZ = static_cast<float>(std::abs(disc.refCoord));
    if (!disc.sensitive || disc.rings.empty()) {
      // Passive, so it has no layers to bound it and carries the radii the
      // reduction found instead.
      subsystem.passives.push_back(PassiveSurfaceDescription{
          .shape = SurfaceShape::Disc,
          .refCoord = absZ,
          .placement = side == SurfaceSide::Positive
                           ? EndcapPlacement::Positive
                           : EndcapPlacement::Negative,
          .minBound = disc.rings.empty() ? static_cast<float>(disc.minBound)
                                         : disc.rings.front().rMin,
          .maxBound = disc.rings.empty() ? static_cast<float>(disc.maxBound)
                                         : disc.rings.back().rMax,
          .material = materialOrNothing(disc)});
      continue;
    }
    endcapFor(subsystem, side)
        .discs.push_back(DiscDescription{.absZ = absZ,
                                         .rings = disc.rings,
                                         .material = materialOrNothing(disc)});
  }

  // Written down rather than left to a position in a list, this being a
  // measurement of a detector: inserting a layer into it later must not re-key
  // the material of the ones behind it.
  assignLayerIndices(description);
  return description;
}

DetectorLayout makeLayoutFromTrackingGeometry(
    const Acts::TrackingGeometry& trackingGeometry,
    const Acts::GeometryContext& gctx,
    const TrackingGeometryLayoutOptions& options) {
  return makeLayout(
      makeDescriptionFromTrackingGeometry(trackingGeometry, gctx, options));
}

std::string defaultSubsystemName(const Acts::Surface& surface) {
  return "volume-" + std::to_string(surface.geometryId().volume());
}

}  // namespace ActsFatras::Synthetic
