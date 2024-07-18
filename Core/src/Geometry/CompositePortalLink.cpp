// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/CompositePortalLink.hpp"

#include <stdexcept>

namespace Acts {

const TrackingVolume* CompositePortalLink::resolveVolume(
    const GeometryContext& /*gctx*/, const Vector2& /*position*/) const {
  throw std::logic_error{"Not implemented"};
}

std::unique_ptr<PortalLinkBase> CompositePortalLink::mergeImpl(
    const PortalLinkBase& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  return other.mergeImpl(*this, surfaceB, surfaceA, direction, logger);
}

std::unique_ptr<PortalLinkBase> CompositePortalLink::mergeImpl(
    const CompositePortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  throw std::logic_error{"Not implemented"};
  ACTS_VERBOSE("Fell through to the Composite mergeImpl with Composite");
}

std::unique_ptr<PortalLinkBase> CompositePortalLink::mergeImpl(
    const TrivialPortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Fell through to the Composite mergeImpl with Trivial");
  throw std::logic_error{"Not implemented"};
}

std::unique_ptr<PortalLinkBase> CompositePortalLink::mergeImpl(
    const GridPortalLink& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Fell through to the Composite mergeImpl with Grid");
  throw std::logic_error{"Not implemented"};
}
}  // namespace Acts
