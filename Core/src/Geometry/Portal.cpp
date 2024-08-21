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
  if (alongNormal == nullptr && oppositeNormal == nullptr) {
    throw std::invalid_argument("At least one link must be provided");
  }

  if (alongNormal != nullptr) {
    setLink(Direction::AlongNormal, std::move(alongNormal));
  }
  if (oppositeNormal != nullptr) {
    setLink(Direction::OppositeNormal, std::move(oppositeNormal));
  }
}

Portal::Portal(Config&& config) {
  if (!config.alongNormal.m_surface && !config.oppositeNormal.m_surface) {
    throw std::invalid_argument("At least one link must be provided");
  }

  if (config.alongNormal.m_surface) {
    setLink(Direction::AlongNormal, std::make_unique<TrivialPortalLink>(
                                        std::move(config.alongNormal.m_surface),
                                        *config.alongNormal.m_volume));
  }
  if (config.oppositeNormal.m_surface) {
    setLink(Direction::OppositeNormal,
            std::make_unique<TrivialPortalLink>(
                std::move(config.oppositeNormal.m_surface),
                *config.oppositeNormal.m_volume));
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
    throw PortalFusingException();
  }
  active = std::move(link);

  // @TODO: To avoid numerical issues with not-exactly-identical surfaces,
  // reset the other side to the exact same surface instance
  if (m_surface == nullptr) {
    m_surface = active->surfacePtr();
  } else {
    // already have a surface, let's set it on the link we just assigned
    active->setSurface(m_surface);
  }
}

void Portal::setLink(Direction direction,
                     std::shared_ptr<RegularSurface> surface,
                     TrackingVolume& volume) {
  setLink(direction,
          std::make_unique<TrivialPortalLink>(std::move(surface), volume));
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

bool Portal::isValid() const {
  return m_alongNormal != nullptr || m_oppositeNormal != nullptr;
}

const RegularSurface& Portal::surface() const {
  assert(m_surface != nullptr);
  return *m_surface;
}

Portal Portal::merge(Portal& aPortal, Portal& bPortal, BinningValue direction,
                     const Logger& logger) {
  ACTS_DEBUG("Merging to portals along " << direction);

  if (&aPortal == &bPortal) {
    ACTS_ERROR("Cannot merge a portal with itself");
    throw PortalMergingException{};
  }

  std::unique_ptr<PortalLinkBase> mergedAlongNormal = nullptr;
  std::unique_ptr<PortalLinkBase> mergedOppositeNormal = nullptr;

  bool aHasAlongNormal = aPortal.m_alongNormal != nullptr;
  bool aHasOppositeNormal = aPortal.m_oppositeNormal != nullptr;
  bool bHasAlongNormal = bPortal.m_alongNormal != nullptr;
  bool bHasOppositeNormal = bPortal.m_oppositeNormal != nullptr;

  if (aHasAlongNormal != bHasAlongNormal ||
      aHasOppositeNormal != bHasOppositeNormal) {
    ACTS_ERROR("Portals do not have the same links attached");
    throw PortalMergingException();
  }

  if (aPortal.m_alongNormal != nullptr) {
    if (bPortal.m_alongNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link along normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links along normal, merging");
    mergedAlongNormal = PortalLinkBase::merge(std::move(aPortal.m_alongNormal),
                                              std::move(bPortal.m_alongNormal),
                                              direction, logger);
  }

  if (aPortal.m_oppositeNormal != nullptr) {
    if (bPortal.m_oppositeNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link opposite normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links opposite normal, merging");
    mergedOppositeNormal = PortalLinkBase::merge(
        std::move(aPortal.m_oppositeNormal),
        std::move(bPortal.m_oppositeNormal), direction, logger);
  }

  aPortal.m_surface.reset();
  bPortal.m_surface.reset();
  return Portal{std::move(mergedAlongNormal), std::move(mergedOppositeNormal)};
}

Portal Portal::fuse(Portal& aPortal, Portal& bPortal, const Logger& logger) {
  ACTS_DEBUG("Fusing two portals");
  if (&aPortal == &bPortal) {
    ACTS_ERROR("Cannot merge a portal with itself");
    throw PortalMergingException{};
  }

  bool aHasAlongNormal = aPortal.m_alongNormal != nullptr;
  bool aHasOppositeNormal = aPortal.m_oppositeNormal != nullptr;
  bool bHasAlongNormal = bPortal.m_alongNormal != nullptr;
  bool bHasOppositeNormal = bPortal.m_oppositeNormal != nullptr;

  if (aPortal.m_surface == nullptr || bPortal.m_surface == nullptr) {
    ACTS_ERROR("Portals have no surface");
    throw PortalFusingException();
  }

  if (*aPortal.m_surface != *bPortal.m_surface) {
    ACTS_ERROR("Portals have different surfaces");
    throw PortalFusingException();
  }

  if (aHasAlongNormal == bHasAlongNormal ||
      aHasOppositeNormal == bHasOppositeNormal) {
    ACTS_ERROR("Portals have the same links attached");
    throw PortalFusingException();
  }

  aPortal.m_surface.reset();
  bPortal.m_surface.reset();
  if (aHasAlongNormal) {
    return Portal{std::move(aPortal.m_alongNormal),
                  std::move(bPortal.m_oppositeNormal)};
  } else {
    return Portal{std::move(bPortal.m_alongNormal),
                  std::move(aPortal.m_oppositeNormal)};
  }
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
    std::unique_ptr<PortalLinkBase> a, std::unique_ptr<PortalLinkBase> b,
    BinningValue direction, const Logger& logger) {
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

  auto gridMerge =
      [&](const GridPortalLink& aGrid,
          const GridPortalLink& bGrid) -> std::unique_ptr<PortalLinkBase> {
    assert(a != nullptr);
    assert(b != nullptr);
    auto merged = GridPortalLink::merge(aGrid, bGrid, direction, logger);
    if (merged == nullptr) {
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);
    }
    return merged;
  };

  if (const auto* aGrid = dynamic_cast<const GridPortalLink*>(a.get());
      aGrid != nullptr) {
    if (const auto* bGrid = dynamic_cast<const GridPortalLink*>(b.get());
        bGrid != nullptr) {
      ACTS_VERBOSE("Merging two grid portals");
      return gridMerge(*aGrid, *bGrid);

    } else if (const auto* bTrivial =
                   dynamic_cast<const TrivialPortalLink*>(b.get());
               bTrivial != nullptr) {
      ACTS_VERBOSE("Merging a grid portal with a trivial portal");
      return gridMerge(*aGrid, *bTrivial->makeGrid(direction));

    } else if (dynamic_cast<const CompositePortalLink*>(b.get()) != nullptr) {
      ACTS_WARNING("Merging a grid portal with a composite portal");
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);

    } else {
      throw std::logic_error{"Portal link type is not supported"};
    }

  } else if (const auto* aTrivial =
                 dynamic_cast<const TrivialPortalLink*>(a.get());
             aTrivial != nullptr) {
    if (const auto* bGrid = dynamic_cast<const GridPortalLink*>(b.get());
        bGrid) {
      ACTS_VERBOSE("Merging a trivial portal with a grid portal");
      return gridMerge(*aTrivial->makeGrid(direction), *bGrid);

    } else if (const auto* bTrivial =
                   dynamic_cast<const TrivialPortalLink*>(b.get());
               bTrivial != nullptr) {
      ACTS_VERBOSE("Merging two trivial portals");
      return gridMerge(*aTrivial->makeGrid(direction),
                       *bTrivial->makeGrid(direction));

    } else if (dynamic_cast<const CompositePortalLink*>(b.get()) != nullptr) {
      ACTS_WARNING("Merging a trivial portal with a composite portal");
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);

    } else {
      throw std::logic_error{"Portal link type is not supported"};
    }

  } else if (dynamic_cast<const CompositePortalLink*>(a.get()) != nullptr) {
    if (dynamic_cast<const GridPortalLink*>(b.get()) != nullptr) {
      ACTS_WARNING("Merging a composite portal with a grid portal");
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);

    } else if (dynamic_cast<const TrivialPortalLink*>(b.get()) != nullptr) {
      ACTS_WARNING("Merging a composite portal with a trivial portal");
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);

    } else if (dynamic_cast<CompositePortalLink*>(b.get()) != nullptr) {
      ACTS_WARNING("Merging two composite portals");
      return std::make_unique<CompositePortalLink>(std::move(a), std::move(b),
                                                   direction);

    } else {
      throw std::logic_error{"Portal link type is not supported"};
    }

  } else {
    throw std::logic_error{"Portal link type is not supported"};
  }
}

}  // namespace Acts
