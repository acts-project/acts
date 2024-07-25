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
    const PortalLinkBase& other, BinningValue direction,
    const Logger& logger) const {
  return other.mergeImpl(*this, direction, logger);
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const CompositePortalLink& other, BinningValue direction,
    const Logger& logger) const {
  (void)other;
  (void)direction;
  (void)logger;
  ACTS_VERBOSE("Merging TrivialPortalLink with CompositePortalLink");
  throw std::logic_error{"Not implemented"};
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const TrivialPortalLink& other, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Merging TrivialPortalLink with TrivialPortalLink");
  ACTS_VERBOSE("Making grids from trivial portal links in " << direction);
  auto grid1 = makeGrid(direction);
  auto grid2 = other.makeGrid(direction);
  return grid1->mergeImpl(*grid2, direction, logger);
}

std::unique_ptr<PortalLinkBase> TrivialPortalLink::mergeImpl(
    const GridPortalLink& other, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Merging TrivialPortalLink with GridPortalLink");
  ACTS_VERBOSE("Making grid from trivial portal link in " << direction);
  auto gridFromTrivial = makeGrid(direction);
  return other.mergeImpl(*gridFromTrivial, direction, logger);
}

}  // namespace Acts
