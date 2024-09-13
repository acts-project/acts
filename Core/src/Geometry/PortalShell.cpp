// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/PortalShell.hpp"

#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <algorithm>
#include <numeric>

namespace Acts {

void CylinderPortalShell::connectOuter(TrackingVolume& volume) {
  for (Face face : {PositiveDisc, NegativeDisc, OuterCylinder, InnerCylinder,
                    NegativePhiPlane, PositivePhiPlane}) {
    auto* portalAtFace = portal(face);
    if (portalAtFace != nullptr) {
      portalAtFace->fill(volume);
    }
  }
}

SingleCylinderPortalShell::SingleCylinderPortalShell(TrackingVolume& volume) {
  if (volume.volumeBounds().type() != VolumeBounds::BoundsType::eCylinder) {
    throw std::invalid_argument("Invalid volume bounds type");
  }

  const auto& bounds =
      dynamic_cast<const CylinderVolumeBounds&>(volume.volumeBounds());

  std::vector<OrientedSurface> orientedSurfaces =
      bounds.orientedSurfaces(volume.transform());

  auto handle = [&](Face face, std::size_t from) {
    const auto& source = orientedSurfaces.at(from);
    m_portals.at(toUnderlying(face)) =
        std::make_shared<Portal>(source.direction, source.surface, volume);
  };

  if (orientedSurfaces.size() == 6) {
    // Fully equipped cylinder
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    handle(InnerCylinder, tubeInnerCover);
    handle(NegativePhiPlane, tubeSectorNegativePhi);
    handle(PositivePhiPlane, tubeSectorPositivePhi);
  } else if (orientedSurfaces.size() == 5) {
    // Phi sector but no inner cylinder (rMin == 0)
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    // Skip inner tube cover, requires offsetting
    handle(NegativePhiPlane, tubeSectorNegativePhi - 1);
    handle(PositivePhiPlane, tubeSectorPositivePhi - 1);
  } else if (orientedSurfaces.size() == 4) {
    // No phi sector but rMin > 0
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    handle(InnerCylinder, tubeInnerCover);
  } else if (orientedSurfaces.size() == 3) {
    // No phi sector and rMin == 0
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
  } else {
    throw std::invalid_argument("Invalid number of oriented surfaces");
  }
}

Portal* SingleCylinderPortalShell::portal(Face face) {
  return portalPtr(face).get();
}

const std::shared_ptr<Portal>& SingleCylinderPortalShell::portalPtr(Face face) {
  return m_portals.at(toUnderlying(face));
}

void SingleCylinderPortalShell::setPortal(std::shared_ptr<Portal> portal,
                                          Face face) {
  m_portals.at(toUnderlying(face)) = std::move(portal);
}

std::size_t SingleCylinderPortalShell::size() const {
  std::size_t count = 0;
  std::ranges::for_each(
      m_portals, [&count](const auto& portal) { count += portal ? 1 : 0; });
  return count;
}

CylinderStackPortalShell::CylinderStackPortalShell(
    const GeometryContext& gctx, std::vector<CylinderPortalShell*> shells,
    BinningValue direction, const Logger& logger)
    : m_direction{direction}, m_shells{std::move(shells)} {
  ACTS_VERBOSE("Making cylinder stack shell in " << m_direction
                                                 << " direction");
  if (std::ranges::any_of(m_shells,
                          [](const auto* shell) { return shell == nullptr; })) {
    ACTS_ERROR("Invalid shell pointer");
    throw std::invalid_argument("Invalid shell pointer");
  }

  auto merge = [&gctx, direction, &shells = m_shells, &logger](Face face) {
    std::vector<std::shared_ptr<Portal>> portals;
    std::transform(shells.begin(), shells.end(), std::back_inserter(portals),
                   [face](auto* shell) { return shell->portalPtr(face); });

    auto merged = std::accumulate(
        std::next(portals.begin()), portals.end(), portals.front(),
        [&](const auto& aPortal,
            const auto& bPortal) -> std::shared_ptr<Portal> {
          assert(aPortal != nullptr);
          assert(bPortal != nullptr);

          return std::make_shared<Portal>(
              Portal::merge(gctx, *aPortal, *bPortal, direction, logger));
        });

    // reset merged portal on all shells
    for (auto& shell : shells) {
      shell->setPortal(merged, face);
    }
  };

  auto fuse = [&gctx, &shells = m_shells, &logger](Face faceA, Face faceB) {
    for (std::size_t i = 1; i < shells.size(); i++) {
      auto& shellA = shells.at(i - 1);
      auto& shellB = shells.at(i);
      auto fused = std::make_shared<Portal>(Portal::fuse(
          gctx, *shellA->portalPtr(faceA), *shellB->portalPtr(faceB), logger));

      shellA->setPortal(fused, faceA);
      shellB->setPortal(fused, faceB);
    }
  };

  if (direction == BinningValue::binR) {
    ACTS_VERBOSE("Merging portals at positive and negative discs");
    merge(PositiveDisc);
    merge(NegativeDisc);

    ACTS_VERBOSE("Fusing portals at outer and inner cylinders");
    fuse(OuterCylinder, InnerCylinder);

  } else if (direction == BinningValue::binZ) {
    bool allHaveInnerCylinders = std::ranges::all_of(
        m_shells, [](const auto* shell) { return shell->size() == 4; });

    bool noneHaveInnerCylinders = std::ranges::all_of(
        m_shells, [](const auto* shell) { return shell->size() == 3; });

    if (!allHaveInnerCylinders && !noneHaveInnerCylinders) {
      ACTS_ERROR("Invalid inner cylinder configuration");
      throw std::invalid_argument("Invalid inner cylinder configuration");
    }

    m_hasInnerCylinder = allHaveInnerCylinders;

    ACTS_VERBOSE("Merging portals at outer cylinders");
    merge(OuterCylinder);

    if (m_hasInnerCylinder) {
      ACTS_VERBOSE("Merging portals at inner cylinders");
      merge(InnerCylinder);
    }

    ACTS_VERBOSE("Fusing portals at positive and negative discs");
    fuse(PositiveDisc, NegativeDisc);

  } else {
    throw std::invalid_argument("Invalid direction");
  }
}

std::size_t CylinderStackPortalShell::size() const {
  return m_hasInnerCylinder ? 4 : 3;
}

Portal* CylinderStackPortalShell::portal(Face face) {
  return portalPtr(face).get();
}

const std::shared_ptr<Portal>& CylinderStackPortalShell::portalPtr(Face face) {
  if (m_direction == BinningValue::binR) {
    switch (face) {
      case NegativeDisc:
        return m_shells.front()->portalPtr(NegativeDisc);
      case PositiveDisc:
        return m_shells.front()->portalPtr(PositiveDisc);
      case OuterCylinder:
        return m_shells.back()->portalPtr(OuterCylinder);
      case InnerCylinder:
        return m_shells.front()->portalPtr(InnerCylinder);
      default:
        std::stringstream ss;
        ss << "Invalid face: " << face;
        throw std::invalid_argument(ss.str());
    }

  } else {
    switch (face) {
      case NegativeDisc:
        return m_shells.front()->portalPtr(NegativeDisc);
      case PositiveDisc:
        return m_shells.back()->portalPtr(PositiveDisc);
      case OuterCylinder:
        [[fallthrough]];
      case InnerCylinder:
        return m_shells.front()->portalPtr(face);
      default:
        std::stringstream ss;
        ss << "Invalid face: " << face;
        throw std::invalid_argument(ss.str());
    }
  }
}

void CylinderStackPortalShell::setPortal(std::shared_ptr<Portal> portal,
                                         Face face) {
  if (m_direction == BinningValue::binR) {
    switch (face) {
      case NegativeDisc:
        [[fallthrough]];
      case PositiveDisc:
        m_shells.front()->setPortal(std::move(portal), face);
        break;
      case OuterCylinder:
        m_shells.back()->setPortal(std::move(portal), OuterCylinder);
        break;
      case InnerCylinder:
        if (!m_hasInnerCylinder) {
          throw std::invalid_argument("Inner cylinder not available");
        }
        m_shells.front()->setPortal(std::move(portal), InnerCylinder);
        break;
      default:
        std::stringstream ss;
        ss << "Invalid face: " << face;
        throw std::invalid_argument(ss.str());
    }

  } else {
    switch (face) {
      case NegativeDisc:
        m_shells.front()->setPortal(std::move(portal), NegativeDisc);
        break;
      case PositiveDisc:
        m_shells.back()->setPortal(std::move(portal), PositiveDisc);
        break;
      case InnerCylinder:
        if (!m_hasInnerCylinder) {
          throw std::invalid_argument("Inner cylinder not available");
        }
        [[fallthrough]];
      case OuterCylinder:
        for (auto* shell : m_shells) {
          shell->setPortal(portal, face);
        }
        break;
      default:
        std::stringstream ss;
        ss << "Invalid face: " << face;
        throw std::invalid_argument(ss.str());
    }
  }
}

std::ostream& operator<<(std::ostream& os, CylinderPortalShell::Face face) {
  switch (face) {
    using enum CylinderPortalShell::Face;
    case PositiveDisc:
      return os << "PositiveDisc";
    case NegativeDisc:
      return os << "NegativeDisc";
    case OuterCylinder:
      return os << "OuterCylinder";
    case InnerCylinder:
      return os << "InnerCylinder";
    case NegativePhiPlane:
      return os << "NegativePhiPlane";
    case PositivePhiPlane:
      return os << "PositivePhiPlane";
    default:
      return os << "Invalid face";
  }
}

}  // namespace Acts
