// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <memory>
#include <stdexcept>

namespace Acts {

Portal::Portal(std::shared_ptr<RegularSurface> surface)
    : m_surface(std::move(surface)) {
  throw_assert(m_surface, "Portal surface is nullptr");
}

const TrackingVolume* Portal::resolveVolume(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction) const {
  const Vector3 normal = m_surface->normal(gctx, position);
  Direction side = Direction::fromScalar(normal.dot(direction));

  const std::unique_ptr<PortalLinkBase>& link =
      side == Direction::AlongNormal ? m_alongNormal : m_oppositeNormal;

  return nullptr;
  if (link == nullptr) {
    // no link is attached in this direction => this is the end of the world as
    // we know it. (i feel fine)
    return nullptr;
  } else {
    // return link->resolveVolume(position);
  }
}

std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link) {
  link.toStream(os);
  return os;
}

void PortalLinkBase::checkMergePreconditions(const PortalLinkBase& a,
                                             const PortalLinkBase& b,
                                             BinningValue direction) {
  const auto& surfaceA = a.surface();
  const auto& surfaceB = b.surface();

  throw_assert(&surfaceA != &surfaceB,
               "Cannot merge portals to the same surface");

  throw_assert(surfaceA.type() == surfaceB.type(),
               "Cannot merge portals of different surface types");

  throw_assert(surfaceA.bounds().type() == surfaceB.bounds().type(),
               "Cannot merge portals of different surface bounds");

  if (const auto* cylA = dynamic_cast<const CylinderSurface*>(&surfaceA);
      cylA != nullptr) {
    const auto* cylB = dynamic_cast<const CylinderSurface*>(&surfaceB);
    throw_assert(cylB != nullptr,
                 "Cannot merge CylinderSurface with "
                 "non-CylinderSurface");
    throw_assert(
        direction == BinningValue::binZ || direction == BinningValue::binRPhi,
        "Invalid binning direction: " + binningValueName(direction));
  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&surfaceA);
             discA != nullptr) {
    const auto* discB = dynamic_cast<const DiscSurface*>(&surfaceB);
    throw_assert(discB != nullptr,
                 "Cannot merge DiscSurface with non-DiscSurface");
    throw_assert(
        direction == BinningValue::binR || direction == BinningValue::binPhi,
        "Invalid binning direction: " + binningValueName(direction));

  } else {
    throw std::logic_error{"Surface type is not supported"};
  }
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
    const std::shared_ptr<const PortalLinkBase>& a,
    const std::shared_ptr<const PortalLinkBase>& b, BinningValue direction,
    const Logger& logger) {
  ACTS_DEBUG("Merging two arbitrary portals");

  ACTS_VERBOSE(" - a:  " << *a);
  ACTS_VERBOSE(" - b: " << *b);

  checkMergePreconditions(*a, *b, direction);

  // Three options:
  // 1. Grid
  // 2. Trivial
  // 3. Composite

  // Grid Grid
  // Grid Trivial
  // Grid Composite
  // Trivial Grid
  // Trivial Trivial
  // Trivial Composite
  // Composite Grid
  // Composite Trivial
  // Composite Composite

  if (auto aGrid = std::dynamic_pointer_cast<const GridPortalLink>(a); aGrid) {
    if (auto bGrid = std::dynamic_pointer_cast<const GridPortalLink>(b);
        bGrid) {
      ACTS_VERBOSE("Merging two grid portals");
      return GridPortalLink::merge(aGrid, bGrid, direction, logger);

    } else if (auto bTrivial =
                   std::dynamic_pointer_cast<const TrivialPortalLink>(b);
               bTrivial) {
      ACTS_WARNING("Merging a grid portal with a trivial portal");
      return GridPortalLink::merge(aGrid, bTrivial->makeGrid(direction),
                                   direction, logger);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<const CompositePortalLink>(b);
               bComposite) {
      ACTS_WARNING("Merging a grid portal with a composite portal");
      return std::make_unique<CompositePortalLink>(aGrid, bComposite,
                                                   direction);

    } else {
      throw std::logic_error{"Portal type is not supported"};
    }

  } else if (auto aTrivial =
                 std::dynamic_pointer_cast<const TrivialPortalLink>(a);
             aTrivial) {
    if (auto bGrid = std::dynamic_pointer_cast<const GridPortalLink>(b);
        bGrid) {
      ACTS_WARNING("Merging a trivial portal with a grid portal");
      return GridPortalLink::merge(aTrivial->makeGrid(direction), bGrid,
                                   direction, logger);

    } else if (auto bTrivial =
                   std::dynamic_pointer_cast<const TrivialPortalLink>(b);
               bTrivial) {
      ACTS_WARNING("Merging two trivial portals");
      return GridPortalLink::merge(aTrivial->makeGrid(direction),
                                   bTrivial->makeGrid(direction), direction,
                                   logger);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<const CompositePortalLink>(b);
               bComposite) {
      ACTS_WARNING("Merging a trivial portal with a composite portal");
      return std::make_unique<CompositePortalLink>(aTrivial, bComposite,
                                                   direction);

    } else {
      throw std::logic_error{"Portal type is not supported"};
    }

  } else if (auto aComposite =
                 std::dynamic_pointer_cast<const CompositePortalLink>(a);
             aComposite) {
    if (auto bGrid = std::dynamic_pointer_cast<const GridPortalLink>(b);
        bGrid) {
      ACTS_WARNING("Merging a composite portal with a grid portal");
      return std::make_unique<CompositePortalLink>(aComposite, bGrid,
                                                   direction);

    } else if (auto bTrivial =
                   std::dynamic_pointer_cast<const TrivialPortalLink>(b);
               bTrivial) {
      ACTS_WARNING("Merging a composite portal with a trivial portal");
      return std::make_unique<CompositePortalLink>(aComposite, bTrivial,
                                                   direction);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<const CompositePortalLink>(b);
               bComposite) {
      ACTS_WARNING("Merging two composite portals");
      return std::make_unique<CompositePortalLink>(aComposite, bComposite,
                                                   direction);

    } else {
      throw std::logic_error{"Portal type is not supported"};
    }

  } else {
    throw std::logic_error{"Portal type is not supported"};
  }
}

}  // namespace Acts
