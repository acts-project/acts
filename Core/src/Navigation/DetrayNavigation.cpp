// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayFwd.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <ranges>
#include <sstream>
#include <stdexcept>

#include <_strings.h>
#include <detray/definitions/grid_axis.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace Acts {

std::unique_ptr<DetraySurfaceGrid> MultiNavigationPolicy::toDetrayPayload(
    const SurfaceLookupFunction& surfaceLookup, const Logger& logger) const {
  // Only ONE of the child policies should return a non-nullptr payload
  for (const auto& policy : m_policyPtrs) {
    auto payload = policy->toDetrayPayload(surfaceLookup, logger);
    if (payload) {
      return payload;
    }
  }
  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid>
Experimental::MultiLayerNavigationPolicy::toDetrayPayload(
    const SurfaceLookupFunction& /*surfaceLookup*/,
    const Logger& /*logger*/) const {
  return nullptr;
}

namespace {

detray::io::axis_payload convertAxis(const IAxis& axis) {
  using enum detray::axis::binning;
  detray::io::axis_payload payload;
  payload.bins = axis.getNBins();
  if (axis.isEquidistant()) {
    payload.binning = e_regular;
    payload.edges = {axis.getMin(), axis.getMax()};
  } else {
    payload.binning = e_irregular;
    payload.edges = axis.getBinEdges();
  }

  switch (axis.getBoundaryType()) {
    using enum Acts::AxisBoundaryType;
    case Open:
      payload.bounds = detray::axis::bounds::e_open;
      break;
    case Closed:
      payload.bounds = detray::axis::bounds::e_circular;
      break;
    case Bound:
      payload.bounds = detray::axis::bounds::e_closed;
      break;
  }

  return payload;
}

detray::axis::label convertAxisDirection(AxisDirection direction) {
  switch (direction) {
    case AxisDirection::AxisX:
      return detray::axis::label::e_x;
    case AxisDirection::AxisY:
      return detray::axis::label::e_y;
    case AxisDirection::AxisZ:
      return detray::axis::label::e_z;
    case AxisDirection::AxisR:
      return detray::axis::label::e_r;
    case AxisDirection::AxisPhi:
      return detray::axis::label::e_phi;
    case AxisDirection::AxisRPhi:
      return detray::axis::label::e_rphi;
    default:
      throw std::invalid_argument(
          "SurfaceArrayNavigationPolicy: Unknown axis direction detected.");
  }
}

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

std::unique_ptr<DetraySurfaceGrid>
SurfaceArrayNavigationPolicy::toDetrayPayload(
    const SurfaceLookupFunction& surfaceLookup, const Logger& logger) const {
  const auto* gridLookup =
      dynamic_cast<const SurfaceArray::ISurfaceGridLookup*>(
          &m_surfaceArray->gridLookup());

  if (gridLookup == nullptr) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "grid based lookup object. This is not currently convertible to "
        "detray");
  }

  AnyGridConstView gridView = [&] {
    auto r = gridLookup->getGridView();
    if (r == std::nullopt) {
      throw std::runtime_error(
          "SurfaceArrayNavigationPolicy: The surface array does not provide a "
          "grid view. This is not currently convertible to detray");
    }
    return r.value();
  }();

  const auto& transform = gridLookup->getTransform();

  constexpr auto tolerance = s_onSurfaceTolerance;

  if ((transform.rotation().matrix() - RotationMatrix3::Identity()).norm() >
      tolerance) {
    throw std::invalid_argument(
        "SurfaceArrayNavigationPolicy: The surface array lookup reports a "
        "rotation. This is not currently convertible to detray");
  }

  ACTS_DEBUG("Converting surface array with " << gridView.dimensions()
                                              << " dims to detray payload");
  std::vector axes = gridLookup->getAxes();

  if (axes.size() != 2) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "2-dimensional grid. This is not currently convertible to detray");
  }

  const IAxis& axis0 = *axes.at(0);
  const IAxis& axis1 = *axes.at(1);

  // Get binning values to determine acceleration structure type
  std::vector binValues = gridLookup->binningValues();
  if (binValues.empty()) {
    // Fall back to default based on surface type
    switch (gridLookup->surfaceType()) {
      using enum Surface::SurfaceType;
      case Cylinder:
        binValues = {AxisDirection::AxisPhi, AxisDirection::AxisZ};
        break;
      case Disc:
        binValues = {AxisDirection::AxisR, AxisDirection::AxisPhi};
        break;
      case Plane:
        binValues = {AxisDirection::AxisX, AxisDirection::AxisY};
        break;
      default:
        throw std::runtime_error(
            "SurfaceArrayNavigationPolicy: Unsupported surface type");
    }
  }

  // Create the detray surface grid payload
  auto gridPayload = std::make_unique<DetraySurfaceGrid>();

  // Set up the grid link with appropriate acceleration structure type
  detray::io::accel_id accelId = getDetrayAccelId(gridLookup->surfaceType());
  gridPayload->grid_link =
      detray::io::typed_link_payload<detray::io::accel_id>{accelId, 0u};

  // @FIXME: We might have to change the order of the axis based on the surface type

  // Convert the axes
  for (std::size_t i = 0; i < axes.size(); ++i) {
    const auto* axis = axes[i];
    ACTS_DEBUG("- Converting axis " << i << " (" << binValues[i]
                                    << "): " << *axis);

    auto axisPayload = convertAxis(*axis);

    // Set axis label based on binning values
    if (i < binValues.size()) {
      axisPayload.label = convertAxisDirection(binValues[i]);
    } else {
      // Default labels if binValues is insufficient
      axisPayload.label =
          (i == 0) ? detray::axis::label::e_x : detray::axis::label::e_y;
    }

    gridPayload->axes.push_back(axisPayload);
  }

  std::set<const Surface*> seenSurfaces;

  using index_type = std::pair<unsigned int, unsigned int>;

  auto fillGeneric = [&](const auto& mapi, const auto& mapj,
                         index_type indices) {
    auto [i, j] = indices;

    auto di = mapi(i);
    auto dj = mapj(j);

    const auto& surfaces = gridView.atLocalBins({i, j});

    std::vector<std::size_t> surfaceIndices;

    for (const auto* surface : surfaces) {
      try {
        std::size_t surfaceIndex = surfaceLookup(surface);
        surfaceIndices.push_back(surfaceIndex);
        seenSurfaces.insert(surface);
      } catch (const std::exception& e) {
        std::stringstream ss;
        ss << "Warning: Could not find surface index for surface "
           << surface->geometryId() << ": " << e.what();
        throw std::runtime_error{ss.str()};
      }
    }

    // Add the bin to the grid payload
    detray::io::grid_bin_payload<std::size_t> binPayload{{{di, dj}},
                                                         surfaceIndices};
    gridPayload->bins.push_back(binPayload);
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

  ACTS_DEBUG("Created detray payload with " << gridPayload->bins.size()
                                            << " populated bins");

  return gridPayload;
}

std::unique_ptr<DetraySurfaceGrid> TryAllNavigationPolicy::toDetrayPayload(
    const SurfaceLookupFunction& /*surfaceLookup*/,
    const Logger& /*logger*/) const {
  return nullptr;
}

}  // namespace Acts
