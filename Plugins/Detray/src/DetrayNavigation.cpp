// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
//
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/DetrayGridExportPolicy.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>

#include <detray/definitions/grid_axis.hpp>
#include <detray/io/frontend/payloads.hpp>

using namespace Acts;

namespace ActsPlugins {

using DetraySurfaceMaterial = DetrayPayloadConverter::DetraySurfaceMaterial;
using DetraySurfaceGrid = DetrayPayloadConverter::DetraySurfaceGrid;

namespace {

detray::io::accel_id getDetrayAccelId(Surface::SurfaceType surfaceType) {
  using enum Surface::SurfaceType;

  switch (surfaceType) {
    case Cylinder:
      return detray::io::accel_id::concentric_cylinder2_grid;
    case Disc:
      return detray::io::accel_id::polar2_grid;
    case Plane:
      return detray::io::accel_id::cartesian2_grid;
    default:
      throw std::runtime_error(
          "SurfaceArrayNavigationPolicy: Unsupported surface type for detray "
          "conversion");
  }
}

}  // namespace

std::optional<DetraySurfaceGrid> DetrayPayloadConverter::convertSurfaceArray(
    const SurfaceArrayNavigationPolicy& policy, const GeometryContext& gctx,
    const SurfaceLookupFunction& surfaceLookup, const Logger& logger) {
  const SurfaceArray& surfaceArray = policy.surfaceArray();

  const Surface* surface = surfaceArray.surfaceRepresentation();
  if (surface == nullptr) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "surface representation. This is not currently convertible to detray");
  }

  const auto& transform = surface->localToGlobalTransform(gctx);

  constexpr auto tolerance = s_onSurfaceTolerance;

  if ((transform.rotation().matrix() - RotationMatrix3::Identity()).norm() >
      tolerance) {
    throw std::invalid_argument(
        "SurfaceArrayNavigationPolicy: The surface array lookup reports a "
        "rotation. This is not currently convertible to detray");
  }

  std::vector axes = surfaceArray.getAxes();
  ACTS_DEBUG("Converting surface array with " << axes.size()
                                              << " dims to detray payload");

  if (axes.size() != 2) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "2-dimensional grid. This is not currently convertible to detray");
  }

  const IAxis& axis0 = *axes.at(0);
  const IAxis& axis1 = *axes.at(1);

  // Get binning values to determine acceleration structure type
  std::vector<AxisDirection> binValues;

  // Fall back to default based on surface type
  switch (surface->type()) {
    using enum Surface::SurfaceType;
    using enum AxisDirection;
    case Cylinder:
      binValues = {AxisPhi, AxisZ};
      break;
    case Disc:
      binValues = {AxisR, AxisPhi};
      break;
    case Plane:
      binValues = {AxisX, AxisY};
      break;
    default:
      throw std::runtime_error(
          "SurfaceArrayNavigationPolicy: Unsupported surface type");
  }

  // Create the detray surface grid payload
  DetraySurfaceGrid gridPayload;

  // Set up the grid link with appropriate acceleration structure type
  detray::io::accel_id accelId = getDetrayAccelId(surface->type());
  gridPayload.grid_link =
      detray::io::typed_link_payload<detray::io::accel_id>{accelId, 0u};

  // @FIXME: We might have to change the order of the axis based on the surface type

  // Convert the axes
  for (std::size_t i = 0; i < axes.size(); ++i) {
    const auto* axis = axes[i];
    ACTS_DEBUG("- Converting axis " << i << " (" << binValues[i]
                                    << "): " << *axis);

    auto axisPayload = DetrayConversionUtils::convertAxis(*axis);

    // Set axis label based on binning values
    if (i < binValues.size()) {
      axisPayload.label =
          DetrayConversionUtils::convertAxisDirection(binValues[i]);
    } else {
      // Default labels if binValues is insufficient
      axisPayload.label =
          (i == 0) ? detray::axis::label::e_x : detray::axis::label::e_y;
    }

    gridPayload.axes.push_back(axisPayload);
  }

  std::set<const Surface*> seenSurfaces;

  using index_type = std::pair<unsigned int, unsigned int>;

  auto fillGeneric = [&](const auto& mapi, const auto& mapj,
                         index_type indices) {
    const auto [i, j] = indices;

    const auto di = mapi(i);
    const auto dj = mapj(j);

    const std::span<const Acts::Surface* const> surfaces =
        surfaceArray.neighbors({i, j}, 0);

    std::vector<std::size_t> surfaceIndices;

    for (const Acts::Surface* srf : surfaces) {
      try {
        std::size_t surfaceIndex = surfaceLookup(srf);
        surfaceIndices.push_back(surfaceIndex);
        seenSurfaces.insert(srf);
      } catch (const std::exception& e) {
        std::stringstream ss;
        ss << "Warning: Could not find surface index for surface "
           << srf->geometryId() << ": " << e.what();
        throw std::runtime_error{ss.str()};
      }
    }

    // Add the bin to the grid payload
    detray::io::grid_bin_payload<std::size_t> binPayload{{{di, dj}},
                                                         surfaceIndices};
    gridPayload.bins.push_back(binPayload);
  };

  // Depending on the axis boundary type, we need to skip over the
  // under/overflow bins, or include them. This needs to be decided by-axis.
  auto makeIndexRange = [](const IAxis& axis) {
    if (axis.getBoundaryType() == AxisBoundaryType::Open) {
      return std::views::iota(0u, axis.getNBins() + 2);
    }
    return std::views::iota(1u, axis.getNBins() + 1);
  };

  auto idx0 = makeIndexRange(axis0);
  auto idx1 = makeIndexRange(axis1);

  auto makeIndexMap =
      [&](const IAxis& axis) -> std::function<unsigned int(unsigned int)> {
    if (axis.getBoundaryType() == AxisBoundaryType::Open) {
      // In case of Open, we loop from [0, N+1], where N is the number of bins.
      // This includes under/overflow bins.
      // Detray also has under/overflow bins in this case, so we keep the
      // indices the same
      return [](unsigned int i) { return i; };
    }

    // For Closed/Bound, detray does not have under/overflow bins.
    // In ACTS, the bins a physically still present, so we only loop
    // [1, N]. Detray's indices however still go [0, N-1], so we subtract 1 from
    // the indices in the direction.
    return [](unsigned int i) { return i - 1; };
  };

  auto fill = std::bind_front(std::bind_front(fillGeneric, makeIndexMap(axis0)),
                              makeIndexMap(axis1));

  for (auto i : idx0) {
    for (auto j : idx1) {
      fill(index_type{i, j});
    }
  }

  ACTS_DEBUG("Filled surfaces " << seenSurfaces.size() << " into grid");

  ACTS_DEBUG("Created detray payload with " << gridPayload.bins.size()
                                            << " populated bins");

  return gridPayload;
}

namespace {

/// Where detray looks up the surface grid of a volume: the virtual reference
/// surface intersected in detray::spatial_grid::search, in the volume frame.
/// Cylinder grids: cylinder at the mean of the smallest and largest portal
/// radius (detray::grid_builder, outer_radius_getter over all portals), which
/// is the outer radius without an inner cylinder portal. Disc grids: plane
/// through the volume centre.
/// @note Must track detray::grid_builder and detray::spatial_grid.
struct DetrayGridReference {
  bool cylinder = true;
  double radius = 0.;
  Transform3 globalToVolume = Transform3::Identity();

  /// Crossing of the line (p, d) with the reference surface, in the volume
  /// frame. Like detray: the nearest solution, unless it lies behind the
  /// track, then the far one; none if both are behind.
  std::optional<Vector3> crossing(const Vector3& globalPosition,
                                  const Vector3& globalDirection) const {
    const Vector3 p = globalToVolume * globalPosition;
    const Vector3 d = globalToVolume.linear() * globalDirection;
    double path = 0.;
    if (cylinder) {
      const double a = d.head<2>().squaredNorm();
      if (a < std::numeric_limits<double>::epsilon()) {
        return std::nullopt;
      }
      const double b = p.head<2>().dot(d.head<2>()) / a;
      const double c = (p.head<2>().squaredNorm() - radius * radius) / a;
      const double discriminant = b * b - c;
      if (discriminant < 0.) {
        return std::nullopt;
      }
      const double s1 = -b + std::sqrt(discriminant);
      const double s2 = -b - std::sqrt(discriminant);
      path = std::abs(s1) < std::abs(s2) ? s1 : s2;
      if (path < 0.) {
        path = std::abs(s1) >= std::abs(s2) ? s1 : s2;
      }
    } else {
      if (std::abs(d.z()) < std::numeric_limits<double>::epsilon()) {
        return std::nullopt;
      }
      path = -p.z() / d.z();
    }
    if (path < 0.) {
      return std::nullopt;
    }
    return p + path * d;
  }
};

std::optional<DetrayGridReference> detrayGridReference(
    const TrackingVolume& volume, const GeometryContext& gctx, bool cylinder) {
  DetrayGridReference ref;
  ref.cylinder = cylinder;
  ref.globalToVolume = volume.localToGlobalTransform(gctx).inverse();
  if (!cylinder) {
    return ref;
  }
  const auto* bounds =
      dynamic_cast<const CylinderVolumeBounds*>(&volume.volumeBounds());
  if (bounds == nullptr) {
    return std::nullopt;
  }
  const double rMin = bounds->get(CylinderVolumeBounds::eMinR);
  const double rMax = bounds->get(CylinderVolumeBounds::eMaxR);
  ref.radius = rMin > 0. ? 0.5 * (rMin + rMax) : rMax;
  return ref;
}

}  // namespace

void DetrayPayloadConverter::fillSurfaceGridForExport(
    DetraySurfaceGrid& grid, const DetrayGridExportPolicy& policy,
    const TrackingVolume& volume, const GeometryContext& gctx,
    const SurfaceLookupFunction& surfaceLookup, const Logger& logger) {
  const DetrayGridExportPolicy::Config& cfg = policy.config();
  const SurfaceArray& surfaceArray = policy.surfaceArrayPolicy().surfaceArray();
  const Surface* representative = surfaceArray.surfaceRepresentation();
  if (representative == nullptr) {
    throw std::runtime_error(
        "DetrayGridExportPolicy: The surface array does not provide a "
        "surface representation");
  }

  const bool isCylinder = representative->type() == Surface::Cylinder;
  if (!isCylinder && representative->type() != Surface::Disc) {
    ACTS_DEBUG(volume.volumeName()
               << ": grid export only fills cylinder and disc grids, keeping "
                  "the plain grid");
    return;
  }

  const std::vector<const IAxis*> axes = surfaceArray.getAxes();
  if (axes.size() != 2) {
    throw std::runtime_error(
        "DetrayGridExportPolicy: The surface array does not provide a "
        "2-dimensional grid");
  }

  const std::optional<DetrayGridReference> ref =
      detrayGridReference(volume, gctx, isCylinder);
  if (!ref.has_value()) {
    ACTS_DEBUG(volume.volumeName()
               << ": cylinder grid in a volume without cylinder volume "
                  "bounds, keeping the plain grid");
    return;
  }

  // The cylinder grid bins phi, the surface measures r * phi
  double phiScale = 1.;
  if (const auto* bounds =
          dynamic_cast<const CylinderBounds*>(&representative->bounds());
      bounds != nullptr) {
    phiScale = bounds->get(CylinderBounds::eR);
  }

  if (isCylinder) {
    ACTS_DEBUG(volume.volumeName()
               << ": layer r = " << phiScale
               << " mm, detray lookup r = " << ref->radius << " mm");
  } else {
    ACTS_DEBUG(volume.volumeName()
               << ": layer z = " << representative->center(gctx).z()
               << " mm, detray lookup at the volume centre");
  }

  using DetrayBin = std::pair<unsigned int, unsigned int>;

  // ACTS bin -> detray bin; same index shift as in convertSurfaceArray
  auto toDetrayBin = [&](std::size_t b0,
                         std::size_t b1) -> std::optional<DetrayBin> {
    const std::array<std::size_t, 2> b{b0, b1};
    std::array<unsigned int, 2> d{};
    for (std::size_t ax = 0; ax < 2; ++ax) {
      if (axes[ax]->getBoundaryType() == AxisBoundaryType::Open) {
        d[ax] = static_cast<unsigned int>(b[ax]);
      } else {
        if (b[ax] < 1 || b[ax] > axes[ax]->getNBins()) {
          return std::nullopt;
        }
        d[ax] = static_cast<unsigned int>(b[ax] - 1);
      }
    }
    return DetrayBin{d[0], d[1]};
  };

  // Detray bin of a point given in the volume frame
  auto detrayBinAt = [&](const Vector3& local) {
    const double phi = VectorHelpers::phi(local);
    return isCylinder
               ? toDetrayBin(axes[0]->getBin(phi), axes[1]->getBin(local.z()))
               : toDetrayBin(axes[0]->getBin(VectorHelpers::perp(local)),
                             axes[1]->getBin(phi));
  };

  // Detray bin -> surfaces the ACTS lookup returns for tracks looked up there
  std::map<DetrayBin, std::set<const Surface*>> lookedUp;
  auto addHits = [&](const std::optional<DetrayBin>& target,
                     std::span<const Surface* const> hits) {
    if (target.has_value()) {
      lookedUp[*target].insert(hits.begin(), hits.end());
    }
  };

  // Fixed evaluation points in grid coordinates: the bin centre and the
  // corners, pulled inside the bin so they cannot fall into a neighbour
  constexpr double kCornerInset = 0.01;
  auto binPoints = [&](std::size_t b0, std::size_t b1) {
    std::array<std::array<double, 2>, 5> points{};
    std::array<std::array<double, 3>, 2> values{};
    for (std::size_t ax = 0; ax < 2; ++ax) {
      const std::size_t bin = ax == 0 ? b0 : b1;
      const double lo = axes[ax]->getBinLowerBound(bin);
      const double hi = axes[ax]->getBinUpperBound(bin);
      const double inset = kCornerInset * (hi - lo);
      values[ax] = {0.5 * (lo + hi), lo + inset, hi - inset};
    }
    points[0] = {values[0][0], values[1][0]};
    points[1] = {values[0][1], values[1][1]};
    points[2] = {values[0][1], values[1][2]};
    points[3] = {values[0][2], values[1][1]};
    points[4] = {values[0][2], values[1][2]};
    return points;
  };

  // Vertex positions bounding the luminous region
  std::vector<double> vertexZs{0.};
  if (cfg.vertexZHalfRange > 0.) {
    vertexZs = {-cfg.vertexZHalfRange, 0., cfg.vertexZHalfRange};
  }

  const SurfaceArray::NeighborWindow window = surfaceArray.neighborWindow();

  for (std::size_t b0 = 1; b0 <= axes[0]->getNBins(); ++b0) {
    for (std::size_t b1 = 1; b1 <= axes[1]->getNBins(); ++b1) {
      if (cfg.fixedWindow) {
        // Direction independent: the widest window ACTS could ever use
        addHits(toDetrayBin(b0, b1),
                surfaceArray.neighbors({b0, b1}, window.max));
        continue;
      }

      for (const auto& [g0, g1] : binPoints(b0, b1)) {
        const Vector2 surfaceLocal =
            isCylinder ? Vector2{g0 * phiScale, g1} : Vector2{g0, g1};
        // The direction only matters for line surfaces
        const Vector3 position =
            representative->localToGlobal(gctx, surfaceLocal, Vector3::UnitZ());

        // Chord-to-tangent angle of the softest track reaching this radius
        double maxBend = 0.;
        if (cfg.minPt > 0.) {
          // In natural units, pT / B is the radius of curvature in mm
          const double curvatureRadius = cfg.minPt / cfg.bField;
          maxBend = std::asin(std::min(
              1., VectorHelpers::perp(position) / (2. * curvatureRadius)));
        }
        std::vector<double> bends{0.};
        if (maxBend > 0.) {
          bends = {-maxBend, 0., maxBend};
        }

        const std::optional<DetrayBin> ownBin =
            detrayBinAt(ref->globalToVolume * position);

        for (const double z0 : vertexZs) {
          const Vector3 straight =
              (position - Vector3{0., 0., z0}).normalized();
          for (const double bend : bends) {
            const Vector3 direction =
                AngleAxis3(bend, Vector3::UnitZ()) * straight;
            // Exactly what ACTS returns, including its dynamic window
            const auto hits = surfaceArray.neighbors(gctx, position, direction);
            if (hits.empty()) {
              continue;
            }
            // Cheap safety: detray falls back to the track position when the
            // reference surface is not ahead
            addHits(ownBin, hits);
            if (const auto crossing = ref->crossing(position, direction);
                crossing.has_value()) {
              addHits(detrayBinAt(*crossing), hits);
            }
          }
        }
      }
    }
  }

  std::size_t entriesBefore = 0;
  std::size_t entriesAfter = 0;
  for (auto& bin : grid.bins) {
    entriesBefore += bin.content.size();
    const auto it = lookedUp.find({bin.loc_index.at(0), bin.loc_index.at(1)});
    if (it != lookedUp.end()) {
      for (const Surface* surface : it->second) {
        bin.content.push_back(surfaceLookup(surface));
      }
      std::ranges::sort(bin.content);
      const auto duplicates = std::ranges::unique(bin.content);
      bin.content.erase(duplicates.begin(), duplicates.end());
    }
    entriesAfter += bin.content.size();
  }

  ACTS_DEBUG(volume.volumeName()
             << ": filled detray grid for export ("
             << (cfg.fixedWindow ? "fixed" : "track") << " window), "
             << grid.bins.size() << " bins, " << entriesBefore << " -> "
             << entriesAfter << " surface entries");
}

#define NOOP_CONVERTER_IMPL_FULL(type, name)                                  \
  std::optional<DetraySurfaceGrid> DetrayPayloadConverter::convert##name(     \
      const type& /*policy*/, const GeometryContext& /*gctx*/,                \
      const SurfaceLookupFunction& /*surfaceLookup*/, const Logger& logger) { \
    ACTS_DEBUG(#name << " does not implement explicit detray conversion");    \
    return std::nullopt;                                                      \
  }

#define NOOP_CONVERTER_IMPL(type) NOOP_CONVERTER_IMPL_FULL(type, type)

NOOP_CONVERTER_IMPL(TryAllNavigationPolicy)
NOOP_CONVERTER_IMPL(MultiNavigationPolicy)
NOOP_CONVERTER_IMPL(CylinderNavigationPolicy)
NOOP_CONVERTER_IMPL(DetrayGridExportPolicy)
NOOP_CONVERTER_IMPL_FULL(MultiLayerNavigationPolicy, MultiLayerNavigationPolicy)

}  // namespace ActsPlugins
