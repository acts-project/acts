// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidPortalShell.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>

#include <boost/algorithm/string/join.hpp>

namespace Acts {

void CuboidPortalShell::fill(TrackingVolume& volume) {
  for (Face face : {negativeXYPlane, positiveXYPlane, negativeYZPlane,
                    positiveYZPlane, negativeZXPlane, positiveZXPlane}) {
    const auto& portalAtFace = portalPtr(face);
    if (portalAtFace != nullptr) {
      portalAtFace->fill(volume);
      volume.addPortal(portalAtFace);
    }
  }
}

SingleCuboidPortalShell::SingleCuboidPortalShell(TrackingVolume& volume)
    : m_volume{&volume} {
  if (m_volume->volumeBounds().type() != VolumeBounds::BoundsType::eCuboid) {
    throw std::invalid_argument("Invalid volume bounds type");
  }

  const auto& bounds =
      dynamic_cast<const CuboidVolumeBounds&>(m_volume->volumeBounds());

  std::vector<OrientedSurface> orientedSurfaces =
      bounds.orientedSurfaces(m_volume->transform());

  auto handle = [&](Face face, std::size_t from) {
    const auto& source = orientedSurfaces.at(from);
    m_portals.at(toUnderlying(face)) =
        std::make_shared<Portal>(source.direction, source.surface, *m_volume);
  };

  handle(negativeXYPlane, negativeFaceXY);
  handle(positiveXYPlane, positiveFaceXY);
  handle(negativeYZPlane, negativeFaceYZ);
  handle(positiveYZPlane, positiveFaceYZ);
  handle(negativeZXPlane, negativeFaceZX);
  handle(positiveZXPlane, positiveFaceZX);
}

Portal* SingleCuboidPortalShell::portal(Face face) {
  return portalPtr(face).get();
}

std::shared_ptr<Portal> SingleCuboidPortalShell::portalPtr(Face face) {
  return m_portals.at(toUnderlying(face));
}

void SingleCuboidPortalShell::setPortal(std::shared_ptr<Portal> portal,
                                        Face face) {
  assert(portal != nullptr);
  assert(portal->isValid());
  m_portals.at(toUnderlying(face)) = std::move(portal);
}

std::size_t SingleCuboidPortalShell::size() const {
  std::size_t count = 0;
  std::ranges::for_each(
      m_portals, [&count](const auto& portal) { count += portal ? 1 : 0; });
  return count;
}

void SingleCuboidPortalShell::applyToVolume() {
  for (std::size_t i = 0; i < m_portals.size(); i++) {
    const auto& portal = m_portals.at(i);
    if (portal != nullptr) {
      if (!portal->isValid()) {
        std::stringstream ss;
        ss << static_cast<Face>(i);
        throw std::runtime_error{"Invalid portal found in shell at " +
                                 ss.str()};
      }
      m_volume->addPortal(portal);
    }
  }
}

bool SingleCuboidPortalShell::isValid() const {
  return std::ranges::all_of(m_portals, [](const auto& portal) {
    return portal == nullptr || portal->isValid();
  });
};

std::string SingleCuboidPortalShell::label() const {
  std::stringstream ss;
  ss << "CuboidShell(vol=" << m_volume->volumeName() << ")";
  return ss.str();
}

CuboidStackPortalShell::CuboidStackPortalShell(
    const GeometryContext& gctx, std::vector<CuboidPortalShell*> shells,
    AxisDirection axis, const Logger& logger)
    : m_axis(axis), m_shells{std::move(shells)} {
  switch (axis) {
    case AxisDirection::AxisX:
      m_direction = Vector3::UnitX();
      break;
    case AxisDirection::AxisY:
      m_direction = Vector3::UnitY();
      break;
    case AxisDirection::AxisZ:
      m_direction = Vector3::UnitZ();
      break;
    default:
      throw std::invalid_argument(axisDirectionName(axis) +
                                  " is not supported ");
  }
  stackShell(gctx, logger);
}

CuboidStackPortalShell::CuboidStackPortalShell(
    const GeometryContext& gctx, std::vector<CuboidPortalShell*> shells,
    const Vector3& direction, const Logger& logger)
    : m_direction{direction}, m_shells{std::move(shells)} {
  Vector3 stackDirection = transform().rotation().inverse() * m_direction;
  m_axis = directionToAxis(stackDirection);
  stackShell(gctx, logger);
}

AxisDirection CuboidStackPortalShell::directionToAxis(
    const Vector3& direction) {
  if ((direction - Vector3::UnitX()).norm() < 1e-4) {
    return AxisDirection::AxisX;
  } else if ((direction - Vector3::UnitY()).norm() < 1e-4) {
    return AxisDirection::AxisY;
  } else if ((direction - Vector3::UnitZ()).norm() < 1e-4) {
    return AxisDirection::AxisZ;
  } else {
    throw std::invalid_argument(
        "CuboidStackPortalShell: Direction does not coincide with axes");
  }
}

void CuboidStackPortalShell::stackShell(const GeometryContext& gctx,
                                        const Logger& logger) {
  ACTS_VERBOSE("Making cuboid stack shell in " << m_axis << " direction");
  if (std::ranges::any_of(m_shells,
                          [](const auto* shell) { return shell == nullptr; })) {
    ACTS_ERROR("Invalid shell pointer");
    throw std::invalid_argument("Invalid shell pointer");
  }

  ACTS_VERBOSE(" ~> " << label());

  if (!std::ranges::all_of(
          m_shells, [](const auto* shell) { return shell->isValid(); })) {
    ACTS_ERROR("Invalid shell");
    throw std::invalid_argument("Invalid shell");
  }
  std::tie(m_frontFace, m_backFace, m_sideFaces) =
      CuboidVolumeBounds::facesFromAxisDirection(m_axis);

  std::unordered_map<Face, AxisDirection> onSurfaceDirs;
  for (Face face : m_sideFaces) {
    const auto& portalAtFace = m_shells.front()->portalPtr(face);
    Vector3 onSurfaceDirection =
        portalAtFace->surface().transform(gctx).rotation().inverse() *
        m_direction;
    onSurfaceDirs[face] = directionToAxis(onSurfaceDirection);
  }

  std::sort(
      m_shells.begin(), m_shells.end(),
      [this](const auto& shellA, const auto& shellB) {
        return (shellA->transform().translation().matrix().dot(m_direction) <
                shellB->transform().translation().matrix().dot(m_direction));
      });

  auto merge = [&](Face face) {
    std::vector<std::shared_ptr<Portal>> portals;
    std::ranges::transform(
        m_shells, std::back_inserter(portals),
        [face](auto* shell) { return shell->portalPtr(face); });

    auto merged = std::accumulate(
        std::next(portals.begin()), portals.end(), portals.front(),
        [&](const auto& aPortal,
            const auto& bPortal) -> std::shared_ptr<Portal> {
          assert(aPortal != nullptr);
          assert(bPortal != nullptr);

          AxisDirection onSurfaceAxis = onSurfaceDirs.at(face);

          return std::make_shared<Portal>(
              Portal::merge(gctx, *aPortal, *bPortal, onSurfaceAxis, logger));
        });

    assert(merged != nullptr);
    assert(merged->isValid());

    // reset merged portal on all shells
    for (auto& shell : m_shells) {
      shell->setPortal(merged, face);
    }
  };

  auto fuse = [&](Face faceA, Face faceB) {
    for (std::size_t i = 1; i < m_shells.size(); i++) {
      auto& shellA = m_shells.at(i - 1);
      auto& shellB = m_shells.at(i);
      ACTS_VERBOSE("Fusing " << shellA->label() << " and " << shellB->label());

      auto fused = std::make_shared<Portal>(Portal::fuse(
          gctx, *shellA->portalPtr(faceA), *shellB->portalPtr(faceB), logger));

      assert(fused != nullptr && "Invalid fused portal");
      assert(fused->isValid() && "Fused portal is invalid");

      shellA->setPortal(fused, faceA);
      shellB->setPortal(fused, faceB);

      assert(shellA->isValid() && "Shell A is not valid after fusing");
      assert(shellB->isValid() && "Shell B is not valid after fusing");
    }
  };

  for (const auto& face : m_sideFaces) {
    merge(face);
  }
  fuse(m_backFace, m_frontFace);
  assert(isValid() && "Shell is not valid after construction");
}

std::size_t CuboidStackPortalShell::size() const {
  return 6;
}

Portal* CuboidStackPortalShell::portal(Face face) {
  return portalPtr(face).get();
}

std::shared_ptr<Portal> CuboidStackPortalShell::portalPtr(Face face) {
  if (face == m_backFace) {
    return m_shells.back()->portalPtr(face);
  } else {
    return m_shells.front()->portalPtr(face);
  }
}

void CuboidStackPortalShell::setPortal(std::shared_ptr<Portal> portal,
                                       Face face) {
  assert(portal != nullptr);

  if (face == m_backFace) {
    m_shells.back()->setPortal(std::move(portal), face);
  } else if (face == m_frontFace) {
    m_shells.front()->setPortal(std::move(portal), face);
  } else {
    for (auto& shell : m_shells) {
      shell->setPortal(portal, face);
    }
  }
}

bool CuboidStackPortalShell::isValid() const {
  return std::ranges::all_of(m_shells, [](const auto* shell) {
    assert(shell != nullptr);
    return shell->isValid();
  });
}

const Transform3& CuboidStackPortalShell::transform() const {
  return m_shells.front()->transform();
}

std::string CuboidStackPortalShell::label() const {
  std::stringstream ss;
  ss << "CuboidStackShell(dir=" << m_axis << ", children=";

  std::vector<std::string> labels;
  std::ranges::transform(m_shells, std::back_inserter(labels),
                         [](const auto* shell) { return shell->label(); });
  ss << boost::algorithm::join(labels, "|");
  ss << ")";
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, CuboidPortalShell::Face face) {
  switch (face) {
    using enum CuboidVolumeBounds::Face;
    case positiveXYPlane:
      return os << "positiveXYPlane";
    case negativeXYPlane:
      return os << "negativeXYPlane";
    case positiveYZPlane:
      return os << "positiveYZPlane";
    case negativeYZPlane:
      return os << "negativeYZPlane";
    case positiveZXPlane:
      return os << "positiveZXPlane";
    case negativeZXPlane:
      return os << "negativeZXPlane";
    default:
      return os << "Invalid face";
  }
}

}  // namespace Acts
