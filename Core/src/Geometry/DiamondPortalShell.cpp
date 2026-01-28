// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DiamondPortalShell.hpp"

#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include <boost/algorithm/string/join.hpp>

namespace Acts {

void DiamondPortalShell::fill(TrackingVolume& volume) {
  for (const auto face : {Face::NegativeZFaceXY, Face::PositiveZFaceXY,
                          Face::NegativeXFaceYZ12, Face::PositiveXFaceYZ12,
                          Face::NegativeXFaceYZ23, Face::PositiveXFaceYZ23,
                          Face::NegativeYFaceZX, Face::PositiveYFaceZX}) {
    const auto& portalAtFace = portalPtr(face);
    if (portalAtFace != nullptr) {
      portalAtFace->fill(volume);
      volume.addPortal(portalAtFace);
    }
  }
}

SingleDiamondPortalShell::SingleDiamondPortalShell(TrackingVolume& volume)
    : m_volume{&volume} {
  if (m_volume->volumeBounds().type() != VolumeBounds::BoundsType::eDiamond) {
    throw std::invalid_argument(
        "SingleDiamondPortalShell: Associated volume does not "
        "have DiamondVolumeBounds");
  }

  const auto& bounds =
      dynamic_cast<const DiamondVolumeBounds&>(m_volume->volumeBounds());

  // fill the protals from the oriented surfaces of the volume bounds
  std::vector<OrientedSurface> orientedSurfaces =
      bounds.boundarySurfaces(*m_volume);
  for (Face face : {Face::NegativeZFaceXY, Face::PositiveZFaceXY,
                    Face::NegativeXFaceYZ12, Face::PositiveXFaceYZ12,
                    Face::NegativeXFaceYZ23, Face::PositiveXFaceYZ23,
                    Face::NegativeYFaceZX, Face::PositiveYFaceZX}) {
    const auto& orientedSurface = orientedSurfaces.at(toUnderlying(face));
    m_portals.at(toUnderlying(face)) = std::make_shared<Portal>(
        orientedSurface.direction, orientedSurface.surface, *m_volume);
  }
}

std::shared_ptr<Portal> SingleDiamondPortalShell::portalPtr(Face face) {
  return m_portals.at(toUnderlying(face));
}

void SingleDiamondPortalShell::setPortal(std::shared_ptr<Portal> portal,
                                         Face face) {
  assert(portal != nullptr);
  assert(portal->isValid());
  m_portals.at(toUnderlying(face)) = std::move(portal);
}

std::size_t SingleDiamondPortalShell::size() const {
  return std::ranges::count_if(
      m_portals, [](const auto& portal) { return portal != nullptr; });
}

void SingleDiamondPortalShell::applyToVolume() {
  for (std::size_t p = 0; p < m_portals.size(); ++p) {
    const auto& portal = m_portals.at(p);
    if (portal != nullptr) {
      if (!portal->isValid()) {
        std::stringstream ss;
        ss << static_cast<Face>(p);
        throw std::runtime_error{"Invalid portal found in shell at" + ss.str()};
      }
      m_volume->addPortal(portal);
    }
  }
}

bool SingleDiamondPortalShell::isValid() const {
  return std::ranges::all_of(m_portals, [](const auto& portal) {
    return portal == nullptr || portal->isValid();
  });
}

std::string SingleDiamondPortalShell::label() const {
  std::stringstream ss;
  ss << "Single Diamond Portal Shell for vol = " + m_volume->volumeName()
     << " ";
  return ss.str();
}

}  // namespace Acts
