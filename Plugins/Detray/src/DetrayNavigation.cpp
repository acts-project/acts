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
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <ranges>
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

/// Where Detray looks a volume's grid up: cylinder grids on the cylinder half
/// way between the volume's portal radii, disc grids on the plane through the
/// volume centre. This is not the ACTS layer surface.
struct DetrayGridReference {
  bool cylinder = true;
  double radius = 0.;
  double planeZ = 0.;
  Transform3 globalToVolume = Transform3::Identity();
};

DetrayGridReference detrayGridReference(const TrackingVolume& volume,
                                        const GeometryContext& gctx,
                                        bool cylinder) {
  DetrayGridReference ref;
  ref.cylinder = cylinder;
  const Transform3& toGlobal = volume.localToGlobalTransform(gctx);
  ref.globalToVolume = toGlobal.inverse();
  ref.planeZ = toGlobal.translation().z();
  if (const auto* bounds =
          dynamic_cast<const CylinderVolumeBounds*>(&volume.volumeBounds());
      bounds != nullptr) {
    const double rMin = bounds->get(CylinderVolumeBounds::eMinR);
    const double rMax = bounds->get(CylinderVolumeBounds::eMaxR);
    // Without an inner portal detray only sees the outer radius
    ref.radius = rMin > 0. ? 0.5 * (rMin + rMax) : rMax;
  }
  return ref;
}

std::vector<double> linspace(double a, double b, unsigned int n) {
  std::vector<double> v;
  if (n <= 1) {
    v.push_back(0.5 * (a + b));
    return v;
  }
  for (unsigned int k = 0; k < n; ++k) {
    v.push_back(a + (b - a) * k / (n - 1));
  }
  return v;
}

}  // namespace

void DetrayPayloadConverter::fillSurfaceGridFromLookup(
    DetraySurfaceGrid& grid, const SurfaceArrayNavigationPolicy& policy,
    const TrackingVolume& volume, const GeometryContext& gctx,
    const SurfaceLookupFunction& surfaceLookup,
    const SurfaceGridFillConfig& cfg, const Logger& logger) {
  const SurfaceArray& surfaceArray = policy.surfaceArray();
  const Surface* representative = surfaceArray.surfaceRepresentation();
  if (representative == nullptr) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "surface representation");
  }

  const bool isCylinder = representative->type() == Surface::Cylinder;
  if (!isCylinder && representative->type() != Surface::Disc) {
    ACTS_DEBUG(volume.volumeName()
               << ": lookup fill only supports cylinder and disc grids, "
                  "keeping the default grid");
    return;
  }

  const std::vector<const IAxis*> axes = surfaceArray.getAxes();
  if (axes.size() != 2) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "2-dimensional grid");
  }

  // A cylinder grid bins phi, but the surface measures r * phi
  double scale0 = 1.;
  if (const auto* cylinder =
          dynamic_cast<const CylinderBounds*>(&representative->bounds());
      cylinder != nullptr) {
    scale0 = cylinder->get(CylinderBounds::eR);
  }

  // Crossing points spread over one ACTS bin, kept off the bin edges
  auto samplePoints = [&](std::size_t ax, std::size_t bin) {
    const double lo = axes[ax]->getBinLowerBound(bin);
    const double hi = axes[ax]->getBinUpperBound(bin);
    const unsigned int n = std::max(cfg.nPointsPerAxis, 1u);
    std::vector<double> points;
    for (unsigned int k = 0; k < n; ++k) {
      points.push_back(lo + (hi - lo) * (k + 0.5) / n);
    }
    return points;
  };
  const std::vector<double> vertexZs =
      linspace(-cfg.vertexZHalfRange, cfg.vertexZHalfRange, cfg.nVertexZ);

  const DetrayGridReference ref = detrayGridReference(volume, gctx, isCylinder);
  if (isCylinder) {
    ACTS_DEBUG(volume.volumeName()
               << ": layer r = " << scale0
               << " mm, detray lookup r = " << ref.radius << " mm");
  } else {
    ACTS_DEBUG(volume.volumeName()
               << ": layer z = " << representative->center(gctx).z()
               << " mm, detray lookup z = " << ref.planeZ << " mm");
  }

  using DetrayBin = std::pair<unsigned int, unsigned int>;
  // ACTS local bins -> detray bin, which has no under/overflow on non-open
  // axes (see convertSurfaceArray)
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
  // Detray bin a global point falls into, in detray's grid frame (the volume)
  auto detrayBinAt = [&](const Vector3& point) {
    const Vector3 local = ref.globalToVolume * point;
    const double phi = std::atan2(local.y(), local.x());
    return ref.cylinder
               ? toDetrayBin(axes[0]->getBin(phi), axes[1]->getBin(local.z()))
               : toDetrayBin(axes[0]->getBin(VectorHelpers::perp(local)),
                             axes[1]->getBin(phi));
  };
  // Crossing of the track line with detray's reference surface, closest to p
  auto referenceCrossing = [&](const Vector3& p,
                               const Vector3& d) -> std::optional<Vector3> {
    if (ref.cylinder) {
      const double a = d.head<2>().squaredNorm();
      const double b = 2. * p.head<2>().dot(d.head<2>());
      const double c = p.head<2>().squaredNorm() - ref.radius * ref.radius;
      const double disc = b * b - 4. * a * c;
      if (a <= 0. || disc < 0.) {
        return std::nullopt;
      }
      const double t1 = (-b - std::sqrt(disc)) / (2. * a);
      const double t2 = (-b + std::sqrt(disc)) / (2. * a);
      return p + (std::abs(t1) < std::abs(t2) ? t1 : t2) * d;
    }
    if (std::abs(d.z()) < 1e-12) {
      return std::nullopt;
    }
    return p + ((ref.planeZ - p.z()) / d.z()) * d;
  };

  // Detray bin -> surfaces that tracks looked up there must see
  std::map<DetrayBin, std::vector<const Surface*>> lookedUp;
  std::vector<std::pair<DetrayBin, const Surface*>> found;
  for (std::size_t b0 = 1; b0 <= axes[0]->getNBins(); ++b0) {
    for (std::size_t b1 = 1; b1 <= axes[1]->getNBins(); ++b1) {
      found.clear();
      for (double g0 : samplePoints(0, b0)) {
        for (double g1 : samplePoints(1, b1)) {
          const Vector3 position = representative->localToGlobal(
              gctx, Vector2{g0 * scale0, g1}, Vector3::UnitZ());
          const double rT = VectorHelpers::perp(position);
          double maxBend = 0.;
          if (cfg.minPt > 0.) {
            // In ACTS natural units, pT / B is the radius of curvature in mm
            const double radiusOfCurvature = cfg.minPt / cfg.bField;
            maxBend = std::asin(std::min(1., rT / (2. * radiusOfCurvature)));
          }
          for (double z0 : vertexZs) {
            const Vector3 straight =
                (position - Vector3{0., 0., z0}).normalized();
            for (double bend : linspace(-maxBend, maxBend, cfg.nBend)) {
              const Vector3 direction =
                  AngleAxis3(bend, Vector3::UnitZ()) * straight;
              const auto hits =
                  surfaceArray.neighbors(gctx, position, direction);
              if (hits.empty()) {
                continue;
              }
              // Detray queries at the reference crossing, or at the track
              // position once the reference is behind it
              std::array<std::optional<DetrayBin>, 2> targets{
                  detrayBinAt(position)};
              if (auto crossing = referenceCrossing(position, direction)) {
                targets[1] = detrayBinAt(*crossing);
              }
              for (const auto& target : targets) {
                if (!target) {
                  continue;
                }
                for (const Surface* surface : hits) {
                  found.emplace_back(*target, surface);
                }
              }
            }
          }
        }
      }
      std::ranges::sort(found);
      const auto last = std::ranges::unique(found);
      found.erase(last.begin(), last.end());
      for (const auto& [target, surface] : found) {
        lookedUp[target].push_back(surface);
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
      const auto last = std::ranges::unique(bin.content);
      bin.content.erase(last.begin(), last.end());
    }
    entriesAfter += bin.content.size();
  }

  ACTS_DEBUG(volume.volumeName()
             << ": " << grid.bins.size() << " bins, " << entriesBefore << " -> "
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
NOOP_CONVERTER_IMPL_FULL(MultiLayerNavigationPolicy, MultiLayerNavigationPolicy)

}  // namespace ActsPlugins
