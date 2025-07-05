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

#include <memory>

#include <detray/definitions/grid_axis.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace Acts {

std::unique_ptr<DetraySurfaceGrid> MultiNavigationPolicy::toDetrayPayload()
    const {
  // Only ONE of the child policies should return a non-nullptr payload
  for (const auto& policy : m_policyPtrs) {
    auto payload = policy->toDetrayPayload();
    if (payload) {
      return payload;
    }
  }
  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid>
Experimental::MultiLayerNavigationPolicy::toDetrayPayload() const {
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

}  // namespace

std::unique_ptr<DetraySurfaceGrid>
SurfaceArrayNavigationPolicy::toDetrayPayload() const {
  const auto* gridLookup =
      dynamic_cast<const SurfaceArray::ISurfaceGridLookup*>(
          &m_surfaceArray->gridLookup());

  if (gridLookup == nullptr) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "grid based lookup object. This is not currently convertible to "
        "detray");
  }

  auto gridView = gridLookup->getGridView();
  if (gridView == std::nullopt) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "grid view. This is not currently convertible to detray");
  }

  const auto& transform = gridLookup->getTransform();

  constexpr auto tolerance = s_onSurfaceTolerance;

  if ((transform.rotation().matrix() - RotationMatrix3::Identity()).norm() >
      tolerance) {
    throw std::invalid_argument(
        "SurfaceArrayNavigationPolicy: The surface array lookup reports a "
        "rotation. This is not currently convertible to detray");
  }

  std::cout << "SurfaceArrayNavigationPolicy: Converting surface array with "
            << gridView->dimensions() << " dims to detray payload" << std::endl;
  auto axes = gridLookup->getAxes();

  if (axes.size() != 2) {
    throw std::runtime_error(
        "SurfaceArrayNavigationPolicy: The surface array does not provide a "
        "2-dimensional grid. This is not currently convertible to detray");
  }

  auto axis1 = convertAxis(*axes[0]);
  auto axis2 = *axes[1];

  DetraySurfaceGrid gridPayload;

  switch (gridLookup->surfaceType()) {
    using enum Surface::SurfaceType;
    case Cylinder:
      gridPayload.axes.emplace_back(convertAxis(const IAxis& axis))
  }

  for (const auto* axis : axes) {
    std::cout << "Axis: " << *axis << std::endl;
    auto& axisPayload = gridPayload.axes.emplace_back(convertAxis(*axis));

    // axisPayload.label = convertAxisDirection();
  }

  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid> TryAllNavigationPolicy::toDetrayPayload()
    const {
  return nullptr;
}

}  // namespace Acts
