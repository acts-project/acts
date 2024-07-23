// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrivialPortalLink.hpp"

#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <memory>

namespace Acts {

std::unique_ptr<GridPortalLink> TrivialPortalLink::makeGrid(
    BinningValue direction) const {
  return GridPortalLink::make(m_surface, *m_volume, direction);
}

const TrackingVolume* TrivialPortalLink::resolveVolume(
    const GeometryContext& /*gctx*/, const Vector2& /*position*/) const {
  return m_volume;
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const PortalLinkBase& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  return other.mergeImpl(*this, surfaceB, surfaceA, direction, logger);
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const CompositePortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  throw std::logic_error{"Not implemented"};
  ACTS_VERBOSE("Fell through to the Trivial mergeImpl with Composite");
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const TrivialPortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Fell through to the Trivial mergeImpl with Trivial");
  throw std::logic_error{"Not implemented"};
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const GridPortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Fell through to the Trivial mergeImpl with Grid");
  throw std::logic_error{"Not implemented"};
}

}  // namespace Acts
