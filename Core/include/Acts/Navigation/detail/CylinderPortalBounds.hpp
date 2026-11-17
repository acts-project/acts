// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace Acts {
class GeometryContext;
class TrackingVolume;
class Portal;
class NavigationStream;
class AppendOnlyNavigationStream;
namespace detail {
/// Immutable supporting faces for the cylindrical portal guard. Stored with the
/// navigation policy, so each propagation only copies the current descriptor.
/// The policy construction context must remain valid for the navigation run.
struct CylinderPortalBounds {
  CylinderPortalBounds() = default;
  CylinderPortalBounds(const GeometryContext& gctx,
                       const TrackingVolume& volume);
  void attach(AppendOnlyNavigationStream& stream) const;
  static const CylinderPortalBounds* get(const NavigationStream& stream) {
    return stream.m_cylinderPortalBounds;
  }
  bool valid() const;
  const TrackingVolume* cachedVolume = nullptr;
  const Portal* lowerDisc = nullptr;
  const Portal* upperDisc = nullptr;
  double lowerZ = 0.;
  double upperZ = 0.;
  const Portal* innerCylinder = nullptr;
  const Portal* outerCylinder = nullptr;
  Vector2 center = Vector2::Zero();
  double innerRadius2 = 0.;
  double outerRadius2 = 0.;
  bool nestedContainer = false;
  std::array<std::pair<const Portal*, double>, 4> containerDiscs{};
  std::size_t containerDiscCount = 0;
  double innerMinZ = 0.;
  double innerMaxZ = 0.;
};
}  // namespace detail
}  // namespace Acts
