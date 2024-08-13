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
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <memory>
#include <stdexcept>

namespace Acts {

Portal::Portal(Direction direction, std::unique_ptr<PortalLinkBase> link) {
  setLink(direction, std::move(link));
  assert(m_surface != nullptr);
}

Portal::Portal(Direction direction, std::shared_ptr<RegularSurface> surface,
               TrackingVolume& volume)
    : Portal(direction,
             std::make_unique<TrivialPortalLink>(std::move(surface), volume)) {}

Portal::Portal(std::unique_ptr<PortalLinkBase> alongNormal,
               std::unique_ptr<PortalLinkBase> oppositeNormal) {
  if (alongNormal != nullptr) {
    setLink(Direction::AlongNormal, std::move(alongNormal));
  }
  if (oppositeNormal != nullptr) {
    setLink(Direction::OppositeNormal, std::move(oppositeNormal));
  }
}

void Portal::setLink(Direction direction,
                     std::unique_ptr<PortalLinkBase> link) {
  assert(link != nullptr);

  auto& active =
      direction == Direction::AlongNormal ? m_alongNormal : m_oppositeNormal;
  auto& other =
      direction == Direction::AlongNormal ? m_oppositeNormal : m_alongNormal;

  // check if surfaces are identical
  if (other != nullptr && link->surface() != other->surface()) {
    throw std::runtime_error("Cannot set two different surfaces");
  }
  active = std::move(link);

  // @TODO: To avoid numerical issues with not-exactly-identical surfaces,
  // reset the other side to the exact same surface instance
  if (m_surface == nullptr) {
    m_surface = active->surfacePtr();
  }
}

const PortalLinkBase* Portal::getLink(Direction direction) const {
  if (direction == Direction::AlongNormal) {
    return m_alongNormal.get();
  } else {
    return m_oppositeNormal.get();
  }
}

const TrackingVolume* Portal::resolveVolume(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction) const {
  assert(m_surface != nullptr);
  const Vector3 normal = m_surface->normal(gctx, position);
  Direction side = Direction::fromScalar(normal.dot(direction));

  const PortalLinkBase* link = side == Direction::AlongNormal
                                   ? m_alongNormal.get()
                                   : m_oppositeNormal.get();

  if (link == nullptr) {
    // no link is attached in this direction => this is the end of the world as
    // we know it. (i feel fine)
    return nullptr;
  } else {
    return link->resolveVolume(gctx, position);
  }
}

std::unique_ptr<Portal> Portal::merge(const std::shared_ptr<Portal>& aPortal,
                                      const std::shared_ptr<Portal>& bPortal,
                                      BinningValue direction,
                                      const Logger& logger) {
  ACTS_DEBUG("Merging to portals along " << direction);
  std::unique_ptr<PortalLinkBase> mergedAlongNormal = nullptr;
  std::unique_ptr<PortalLinkBase> mergedOppositeNormal = nullptr;

  bool aHasAlongNormal = aPortal->m_alongNormal != nullptr;
  bool aHasOppositeNormal = aPortal->m_oppositeNormal != nullptr;
  bool bHasAlongNormal = bPortal->m_alongNormal != nullptr;
  bool bHasOppositeNormal = bPortal->m_oppositeNormal != nullptr;

  if (aHasAlongNormal != bHasAlongNormal ||
      aHasOppositeNormal != bHasOppositeNormal) {
    ACTS_ERROR("Portals do not have the same links attached");
    throw PortalMergingException();
  }

  if (aPortal->m_alongNormal != nullptr) {
    if (bPortal->m_alongNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link along normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links along normal, merging");
    mergedAlongNormal = PortalLinkBase::merge(
        aPortal->m_alongNormal, bPortal->m_alongNormal, direction, logger);
  }

  if (aPortal->m_oppositeNormal != nullptr) {
    if (bPortal->m_oppositeNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link opposite normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links opposite normal, merging");
    mergedOppositeNormal =
        PortalLinkBase::merge(aPortal->m_oppositeNormal,
                              bPortal->m_oppositeNormal, direction, logger);
  }

  return std::make_unique<Portal>(std::move(mergedAlongNormal),
                                  std::move(mergedOppositeNormal));
}

std::unique_ptr<Portal> Portal::fuse(const std::shared_ptr<Portal>& aPortal,
                                     const std::shared_ptr<Portal>& bPortal,
                                     const Logger& logger) {
  return nullptr;
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

    throw_assert(dynamic_cast<const RadialBounds*>(&discA->bounds()) &&
                     dynamic_cast<const RadialBounds*>(&discB->bounds()),
                 "DiscSurface bounds must be RadialBounds");

  } else {
    throw std::logic_error{"Surface type is not supported"};
  }
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
    const std::shared_ptr<PortalLinkBase>& a,
    const std::shared_ptr<PortalLinkBase>& b, BinningValue direction,
    const Logger& logger) {
  ACTS_DEBUG("Merging two arbitrary portals");

  ACTS_VERBOSE(" - a: " << *a);
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

  if (auto aGrid = std::dynamic_pointer_cast<GridPortalLink>(a); aGrid) {
    if (auto bGrid = std::dynamic_pointer_cast<GridPortalLink>(b); bGrid) {
      ACTS_VERBOSE("Merging two grid portals");
      return GridPortalLink::merge(aGrid, bGrid, direction, logger);

    } else if (auto bTrivial = std::dynamic_pointer_cast<TrivialPortalLink>(b);
               bTrivial) {
      ACTS_VERBOSE("Merging a grid portal with a trivial portal");
      return GridPortalLink::merge(aGrid, bTrivial->makeGrid(direction),
                                   direction, logger);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<CompositePortalLink>(b);
               bComposite) {
      ACTS_WARNING("Merging a grid portal with a composite portal");
      return std::make_unique<CompositePortalLink>(aGrid, bComposite,
                                                   direction);

    } else {
      throw std::logic_error{"Portal type is not supported"};
    }

  } else if (auto aTrivial = std::dynamic_pointer_cast<TrivialPortalLink>(a);
             aTrivial) {
    if (auto bGrid = std::dynamic_pointer_cast<GridPortalLink>(b); bGrid) {
      ACTS_VERBOSE("Merging a trivial portal with a grid portal");
      return GridPortalLink::merge(aTrivial->makeGrid(direction), bGrid,
                                   direction, logger);

    } else if (auto bTrivial =
                   std::dynamic_pointer_cast<const TrivialPortalLink>(b);
               bTrivial) {
      ACTS_VERBOSE("Merging two trivial portals");
      return GridPortalLink::merge(aTrivial->makeGrid(direction),
                                   bTrivial->makeGrid(direction), direction,
                                   logger);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<CompositePortalLink>(b);
               bComposite) {
      ACTS_WARNING("Merging a trivial portal with a composite portal");
      return std::make_unique<CompositePortalLink>(aTrivial, bComposite,
                                                   direction);

    } else {
      throw std::logic_error{"Portal type is not supported"};
    }

  } else if (auto aComposite =
                 std::dynamic_pointer_cast<CompositePortalLink>(a);
             aComposite) {
    if (auto bGrid = std::dynamic_pointer_cast<GridPortalLink>(b); bGrid) {
      ACTS_WARNING("Merging a composite portal with a grid portal");
      return std::make_unique<CompositePortalLink>(aComposite, bGrid,
                                                   direction);

    } else if (auto bTrivial = std::dynamic_pointer_cast<TrivialPortalLink>(b);
               bTrivial) {
      ACTS_WARNING("Merging a composite portal with a trivial portal");
      return std::make_unique<CompositePortalLink>(aComposite, bTrivial,
                                                   direction);

    } else if (auto bComposite =
                   std::dynamic_pointer_cast<CompositePortalLink>(b);
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
