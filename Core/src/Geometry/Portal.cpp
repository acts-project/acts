// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
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

}  // namespace Acts
