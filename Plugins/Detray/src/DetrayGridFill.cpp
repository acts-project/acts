// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <detray/definitions/grid_axis.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace ActsPlugins {

namespace {

using Acts::Vector3;

/// One axis of the detray grid, binned the way detray bins it
struct DetrayAxis {
  double min = 0.;
  double width = 1.;
  int nBins = 0;
  bool circular = false;

  /// Unwrapped bin index, as detray's regular binning computes it
  int rawBin(double v) const {
    return static_cast<int>((v - min) / width + 1.) - 1;
  }

  /// Bin index after the boundary treatment of the axis
  int bin(int raw) const {
    if (circular) {
      return (nBins + (raw % nBins)) % nBins;
    }
    return std::clamp(raw, 0, nBins - 1);
  }
};

/// Helix from a vertex on the beam axis, parametrised by transverse arc
/// length s; a curvature of zero is a straight line
struct Helix {
  double z0 = 0.;
  double kappa = 0.;
  double psi0 = 0.;
  double cotTheta = 0.;

  Vector3 position(double s) const {
    const double half = 0.5 * kappa * s;
    const double chord =
        std::abs(half) < 1e-9 ? s : 2. * std::sin(half) / kappa;
    return {chord * std::cos(psi0 + half), chord * std::sin(psi0 + half),
            z0 + cotTheta * s};
  }

  Vector3 direction(double s) const {
    const double psi = psi0 + kappa * s;
    return Vector3{std::cos(psi), std::sin(psi), cotTheta}.normalized();
  }

  /// Arc length at which the transverse distance to the axis is @p r
  std::optional<double> arcAtRadius(double r) const {
    if (std::abs(kappa) < 1e-12) {
      return r;
    }
    const double arg = 0.5 * std::abs(kappa) * r;
    if (arg > 1.) {
      return std::nullopt;
    }
    return 2. * std::asin(arg) / std::abs(kappa);
  }

  /// Helix through the vertex (0, 0, z0) and @p point with curvature @p k
  static std::optional<std::pair<Helix, double>> through(const Vector3& point,
                                                         double z0, double k) {
    Helix h;
    h.z0 = z0;
    h.kappa = k;
    const double rho = std::hypot(point.x(), point.y());
    const auto sPoint = h.arcAtRadius(rho);
    if (!sPoint.has_value() || *sPoint <= 0.) {
      return std::nullopt;
    }
    h.cotTheta = (point.z() - z0) / *sPoint;
    h.psi0 = std::atan2(point.y(), point.x()) - 0.5 * k * *sPoint;
    return std::pair{h, *sPoint};
  }
};

}  // namespace

std::optional<DetrayPayloadConverter::DetraySurfaceGrid>
DetrayPayloadConverter::convertSurfaceGridForDetray(
    const Acts::SurfaceArrayNavigationPolicy& policy,
    const Acts::TrackingVolume& volume, const Acts::GeometryContext& gctx,
    const SurfaceLookupFunction& surfaceLookup,
    const DetraySurfaceGridConfig& fill, const Acts::Logger& logger) {
  using Acts::SurfaceArrayNavigationPolicy;

  // Rebin the axes the surface array bins automatically; detray pays for every
  // copy of a surface, so it wants coarser bins than ACTS navigation
  const SurfaceArrayNavigationPolicy::Config& actsConfig = policy.config();
  std::optional<SurfaceArrayNavigationPolicy> rebinned;
  if (fill.numberOfBinsFactor.has_value() &&
      *fill.numberOfBinsFactor != actsConfig.numberOfBinsFactor &&
      (actsConfig.bins.first == 0 || actsConfig.bins.second == 0)) {
    SurfaceArrayNavigationPolicy::Config config = actsConfig;
    config.numberOfBinsFactor = *fill.numberOfBinsFactor;
    rebinned.emplace(gctx, volume, logger, config);
  }
  const SurfaceArrayNavigationPolicy& binning =
      rebinned.has_value() ? *rebinned : policy;

  std::optional<DetraySurfaceGrid> grid =
      convertSurfaceArray(binning, gctx, surfaceLookup, logger);
  if (grid.has_value()) {
    fillSurfaceGridDetrayNative(*grid, binning.surfaceArray(), volume, gctx,
                                surfaceLookup, fill, logger);
  }
  return grid;
}

void DetrayPayloadConverter::fillSurfaceGridDetrayNative(
    DetraySurfaceGrid& grid, const Acts::SurfaceArray& surfaceArray,
    const Acts::TrackingVolume& volume, const Acts::GeometryContext& gctx,
    const SurfaceLookupFunction& surfaceLookup,
    const DetraySurfaceGridConfig& fill, const Acts::Logger& logger) {
  using namespace Acts;

  const Surface* representative = surfaceArray.surfaceRepresentation();
  const bool isCylinder =
      representative != nullptr && representative->type() == Surface::Cylinder;
  const bool isDisc =
      representative != nullptr && representative->type() == Surface::Disc;
  const auto* bounds =
      dynamic_cast<const CylinderVolumeBounds*>(&volume.volumeBounds());
  if ((!isCylinder && !isDisc) || bounds == nullptr || grid.axes.size() != 2) {
    ACTS_WARNING(volume.volumeName()
                 << ": no detray-native fill for this grid, keeping it");
    return;
  }

  const Transform3& volumeTransform = volume.localToGlobalTransform(gctx);
  const Vector3 volumeCenter = volumeTransform.translation();
  if ((volumeTransform.rotation().matrix() - RotationMatrix3::Identity())
              .norm() > 1e-9 ||
      std::hypot(volumeCenter.x(), volumeCenter.y()) > 1e-6) {
    ACTS_WARNING(volume.volumeName()
                 << ": volume is rotated or off the beam axis, keeping the "
                    "plain grid");
    return;
  }

  const double rMin = bounds->get(CylinderVolumeBounds::eMinR);
  const double rMax = bounds->get(CylinderVolumeBounds::eMaxR);
  const double halfZ = bounds->get(CylinderVolumeBounds::eHalfLengthZ);
  const double zLow = volumeCenter.z() - halfZ;
  const double zHigh = volumeCenter.z() + halfZ;
  // Mirrors detray::grid_builder: mean of the smallest and largest portal
  // radius, i.e. the outer radius without an inner cylinder portal
  double refRadius = rMin > 0. ? 0.5 * (rMin + rMax) : rMax;
  if (isCylinder && fill.lookupOnLayer && !surfaceArray.surfaces().empty()) {
    double sum = 0.;
    for (const Surface* surface : surfaceArray.surfaces()) {
      const Vector3 r = surface->referencePosition(gctx, AxisDirection::AxisR);
      sum += std::hypot(r.x(), r.y());
    }
    refRadius = sum / static_cast<double>(surfaceArray.surfaces().size());
    grid.reference_radius = static_cast<detray::io::scalar>(refRadius);
  }

  std::array<DetrayAxis, 2> axes{};
  for (std::size_t i = 0; i < 2; ++i) {
    const auto& ax = grid.axes[i];
    if (ax.binning != detray::axis::binning::e_regular ||
        ax.edges.size() != 2 || ax.bins == 0) {
      ACTS_WARNING(volume.volumeName()
                   << ": irregular grid axis, keeping the plain grid");
      return;
    }
    axes[i].min = ax.edges[0];
    axes[i].nBins = static_cast<int>(ax.bins);
    axes[i].width = (ax.edges[1] - ax.edges[0]) / static_cast<double>(ax.bins);
    axes[i].circular = ax.bounds == detray::axis::bounds::e_circular;
  }
  // Cylinder grids bin (phi, z), disc grids (r, phi)
  const std::size_t phiAxis = isCylinder ? 0 : 1;
  const std::size_t spanAxis = 1 - phiAxis;

  // detray's lookup point in grid coordinates: the tangent at the query
  // position crosses the virtual reference surface; detray takes the nearest
  // crossing even if it lies behind (overstep tolerance is -max)
  auto lookupPoint = [&](const Vector3& pos,
                         const Vector3& dir) -> std::array<double, 2> {
    if (isCylinder) {
      const double dPerp2 = dir.x() * dir.x() + dir.y() * dir.y();
      const double pPerp2 = pos.x() * pos.x() + pos.y() * pos.y();
      const double radDiff = refRadius * refRadius - pPerp2;
      double path = 0.;
      if (dPerp2 < std::numeric_limits<double>::epsilon()) {
        return {std::atan2(pos.y(), pos.x()), pos.z()};
      }
      if (std::abs(radDiff) > 1. * UnitConstants::um) {
        const double k = -(pos.x() * dir.x() + pos.y() * dir.y()) / dPerp2;
        const double discr = radDiff / dPerp2 + k * k;
        if (discr < 0.) {
          return {std::atan2(pos.y(), pos.x()), pos.z()};
        }
        const double s1 = k + std::sqrt(discr);
        const double s2 = k - std::sqrt(discr);
        path = std::abs(s1) < std::abs(s2) ? s1 : s2;
      }
      const Vector3 x = pos + path * dir;
      return {std::atan2(x.y(), x.x()), x.z()};
    }
    // Disc grid: plane through the volume centre
    Vector3 x = pos;
    if (std::abs(dir.z()) > std::numeric_limits<double>::epsilon()) {
      x = pos + ((volumeCenter.z() - pos.z()) / dir.z()) * dir;
    }
    return {std::hypot(x.x(), x.y()), std::atan2(x.y(), x.x())};
  };

  // Track family to serve
  std::vector<double> vertexZs;
  const std::size_t nZ = std::max<std::size_t>(fill.nVertexZ, 1);
  for (std::size_t i = 0; i < nZ; ++i) {
    vertexZs.push_back(nZ == 1 ? 0.
                               : fill.vertexZHalfRange *
                                     (-1. + 2. * static_cast<double>(i) /
                                                static_cast<double>(nZ - 1)));
  }
  std::vector<double> curvatures{0.};
  if (fill.minPt > 0. && fill.bField > 0. && fill.nCurvature > 1) {
    // Radius of curvature in mm: R = pT / (0.3 B), pT in GeV and B in T
    const double radius = (fill.minPt / UnitConstants::GeV) /
                          (0.299792458 * fill.bField / UnitConstants::T) *
                          UnitConstants::m;
    const double kMax = 1. / radius;
    curvatures.clear();
    const std::size_t nK = fill.nCurvature;
    for (std::size_t i = 0; i < nK; ++i) {
      curvatures.push_back(kMax * (-1. + 2. * static_cast<double>(i) /
                                             static_cast<double>(nK - 1)));
    }
  }

  // Query position of a track hitting @p point: where it enters the volume
  auto entry =
      [&](const Helix& h,
          double sPoint) -> std::optional<std::pair<Vector3, Vector3>> {
    double sEnter = 0.;
    if (rMin > 0.) {
      const auto s = h.arcAtRadius(rMin);
      if (!s.has_value()) {
        return std::nullopt;
      }
      sEnter = std::max(sEnter, *s);
    }
    if (h.z0 < zLow || h.z0 > zHigh) {
      if (std::abs(h.cotTheta) < 1e-12) {
        return std::nullopt;
      }
      const double zFace = h.z0 < zLow ? zLow : zHigh;
      const double s = (zFace - h.z0) / h.cotTheta;
      if (s < 0.) {
        return std::nullopt;
      }
      sEnter = std::max(sEnter, s);
    }
    if (sEnter > sPoint * (1. + 1e-9) + 1e-6) {
      return std::nullopt;
    }
    return std::pair{h.position(sEnter), h.direction(sEnter)};
  };

  // detray bin -> surfaces
  std::map<std::pair<int, int>, std::set<std::size_t>> contents;
  std::size_t nUnreached = 0;

  const std::vector<const Surface*> surfaces = surfaceArray.surfaces();
  for (const Surface* surface : surfaces) {
    const std::size_t index = surfaceLookup(surface);
    // Closed outline of the surface, densely sampled below
    const Polyhedron polyhedron = surface->polyhedronRepresentation(gctx, 16);

    // Unwrap phi around the outline: the surface centre can sit far from the
    // sensor, e.g. the origin of an annulus
    Vector3 centroid = Vector3::Zero();
    for (const Vector3& v : polyhedron.vertices) {
      centroid += v;
    }
    const double centerPhi = std::atan2(centroid.y(), centroid.x());

    bool reached = false;
    for (const double z0 : vertexZs) {
      for (const double kappa : curvatures) {
        // Per column of the span axis, the phi range of the outline,
        // unwrapped around the surface
        std::map<int, std::pair<double, double>> columns;
        auto addSample = [&](const Vector3& point) {
          const auto helix = Helix::through(point, z0, kappa);
          if (!helix.has_value()) {
            return std::optional<std::array<double, 2>>{};
          }
          const auto query = entry(helix->first, helix->second);
          if (!query.has_value()) {
            return std::optional<std::array<double, 2>>{};
          }
          std::array<double, 2> g = lookupPoint(query->first, query->second);
          double dphi = g[phiAxis] - centerPhi;
          dphi = std::remainder(dphi, 2. * std::numbers::pi);
          g[phiAxis] = centerPhi + dphi;
          const int col =
              axes[spanAxis].bin(axes[spanAxis].rawBin(g[spanAxis]));
          const double phi = g[phiAxis];
          auto [it, inserted] = columns.try_emplace(col, phi, phi);
          if (!inserted) {
            it->second.first = std::min(it->second.first, phi);
            it->second.second = std::max(it->second.second, phi);
          }
          return std::optional<std::array<double, 2>>{g};
        };

        for (const auto& face : polyhedron.faces) {
          for (std::size_t i = 0; i < face.size(); ++i) {
            const Vector3& from = polyhedron.vertices.at(face.at(i));
            const Vector3& to =
                polyhedron.vertices.at(face.at((i + 1) % face.size()));
            // Subdivide until neighbouring samples are close in bins
            const auto a = addSample(from);
            const auto b = addSample(to);
            if (!a.has_value() || !b.has_value()) {
              continue;
            }
            const double spanBins =
                std::max(std::abs((*b)[0] - (*a)[0]) / axes[0].width,
                         std::abs((*b)[1] - (*a)[1]) / axes[1].width);
            const auto nSub = static_cast<std::size_t>(
                std::ceil(spanBins / fill.maxSampleSpacing));
            for (std::size_t k = 1; k < nSub; ++k) {
              const double t =
                  static_cast<double>(k) / static_cast<double>(nSub);
              addSample(from + t * (to - from));
            }
          }
        }

        // detray bins the atan2 phi in [-pi, pi]. The phi axis need not
        // span 2 pi, so the parts of a range beyond the seam are binned from
        // their wrapped values rather than by wrapping bin indices.
        constexpr double pi = std::numbers::pi;
        auto fillRows = [&](int col, double lo, double hi) {
          if (lo > hi) {
            return;
          }
          for (int row = axes[phiAxis].rawBin(lo);
               row <= axes[phiAxis].rawBin(hi); ++row) {
            const int rowBin = axes[phiAxis].bin(row);
            const std::pair<int, int> key =
                phiAxis == 0 ? std::pair{rowBin, col} : std::pair{col, rowBin};
            contents[key].insert(index);
            reached = true;
          }
        };
        for (const auto& [col, range] : columns) {
          const auto [lo, hi] = range;
          fillRows(col, std::max(lo, -pi), std::min(hi, pi));
          if (hi > pi) {
            fillRows(col, -pi, hi - 2. * pi);
          }
          if (lo < -pi) {
            fillRows(col, lo + 2. * pi, pi);
          }
        }
      }
    }
    if (!reached) {
      ++nUnreached;
    }
  }

  std::size_t entriesBefore = 0;
  std::size_t entriesAfter = 0;
  std::size_t maxContent = 0;
  for (auto& bin : grid.bins) {
    entriesBefore += bin.content.size();
    const auto it = contents.find({static_cast<int>(bin.loc_index.at(0)),
                                   static_cast<int>(bin.loc_index.at(1))});
    // Keep the ACTS assignment where no track of the family looks
    if (it != contents.end()) {
      bin.content.assign(it->second.begin(), it->second.end());
    }
    entriesAfter += bin.content.size();
    maxContent = std::max(maxContent, bin.content.size());
  }

  if (nUnreached > 0) {
    ACTS_WARNING(volume.volumeName()
                 << ": " << nUnreached
                 << " surfaces are not reached by the track family");
  }
  ACTS_DEBUG(volume.volumeName()
             << ": detray-native fill, reference r = " << refRadius << " mm, "
             << grid.bins.size() << " bins, entries " << entriesBefore << " -> "
             << entriesAfter << ", max per bin " << maxContent);
}

}  // namespace ActsPlugins
