// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/TrackingGeometryPrintVisitor.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

namespace {
inline std::string whiteSpaces(const std::size_t n) {
  std::string str{};
  str.resize(n, ' ');
  return str;
}
}  // namespace

namespace Acts::detail {
TrackingGeometryPrintVisitor::TrackingGeometryPrintVisitor(
    const Acts::GeometryContext& gctx, std::size_t indentation)
    : m_gctx{gctx}, m_indentation{indentation} {}
TrackingGeometryPrintVisitor::~TrackingGeometryPrintVisitor() = default;
std::stringstream& TrackingGeometryPrintVisitor::stream() {
  return m_printStream;
}
const std::stringstream& TrackingGeometryPrintVisitor::stream() const {
  return m_printStream;
}
void TrackingGeometryPrintVisitor::visitVolume(
    const Acts::TrackingVolume& volume) {
  updateDepth(volume);
  m_printStream << whiteSpaces(m_currentDepth) << volNumber(volume) << ") "
                << (volume.isAlignable() ? "Alignable volume " : "Volume ")
                << volume.volumeName() << " @ "
                << toString(volume.center(m_gctx))
                << " --- id: " << volume.geometryId()
                << " #surfaces: " << volume.surfaces().size()
                << ", #portals: " << volume.portals().size()
                << ", #sub-volumes: " << volume.volumes().size() << std::endl;
}
std::size_t TrackingGeometryPrintVisitor::volNumber(
    const TrackingVolume& trkVol) const {
  const TrackingVolume* parent = trkVol.motherVolume();
  if (parent == nullptr) {
    return 1;
  }
  return 1 + std::distance(
                 parent->volumes().begin(),
                 std::ranges::find_if(parent->volumes(),
                                      [&trkVol](const TrackingVolume& child) {
                                        return &child == &trkVol;
                                      }));
}
void TrackingGeometryPrintVisitor::visitPortal(const Portal& portal) {
  const auto& surf = portal.surface();
  m_printStream << whiteSpaces(m_currentDepth + m_indentation) << " ++++ "
                << Surface::s_surfaceTypeNames[toUnderlying(surf.type())] <<

      " portal  --- id: " << surf.geometryId() << ", bounds: " << surf.bounds()
                << ", alignable: " << (surf.isAlignable() ? "yay" : "nay")
                << ", material: "
                << (surf.surfaceMaterial() != nullptr ? "yay" : "nay")
                << std::endl;
}

void TrackingGeometryPrintVisitor::visitSurface(const Surface& surface) {
  m_printStream << whiteSpaces(m_currentDepth + m_indentation) << " **** "
                << Surface::s_surfaceTypeNames[toUnderlying(surface.type())]
                << " surface "
                << " @ " << toString(surface.center(m_gctx))
                << " --- id: " << surface.geometryId()
                << ", sensitive: " << (surface.isSensitive() ? "yay" : "nay")
                << ", alignable: " << (surface.isAlignable() ? "yay" : "nay")
                << ", material: "
                << (surface.surfaceMaterial() != nullptr ? "yay" : "nay")
                << std::endl;
}

void TrackingGeometryPrintVisitor::updateDepth(const TrackingVolume& trkVol) {
  m_currentDepth = 0;
  const TrackingVolume* parent = trkVol.motherVolume();
  while (parent != nullptr) {
    ++m_currentDepth;
    parent = parent->motherVolume();
  }
  m_currentDepth *= m_indentation;
}

}  // namespace Acts::detail
