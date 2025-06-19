// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"

#include "Acts/Surfaces/Surface.hpp"

#include <utility>

namespace Acts {

Geant4DetectorElement::Geant4DetectorElement(std::shared_ptr<Surface> surface,
                                             const G4VPhysicalVolume& g4physVol,
                                             const Transform3& toGlobal,
                                             double thickness)
    : m_surface(std::move(surface)),
      m_g4physVol(&g4physVol),
      m_toGlobal(toGlobal),
      m_thickness(thickness) {
  if (m_surface == nullptr) {
    throw std::invalid_argument(
        "Geant4DetectorElement: Surface cannot be nullptr");
  }
  if (m_surface->associatedDetectorElement() != nullptr) {
    throw std::logic_error(
        "Geant4DetectorElement: Surface already has an associated detector "
        "element");
  }
  m_surface->assignDetectorElement(*this);
}

const Transform3& Geant4DetectorElement::transform(
    const GeometryContext& /*gctx*/) const {
  return m_toGlobal;
}

const Surface& Geant4DetectorElement::surface() const {
  return *m_surface;
}

Surface& Geant4DetectorElement::surface() {
  return *m_surface;
}

double Geant4DetectorElement::thickness() const {
  return m_thickness;
}

const G4VPhysicalVolume& Geant4DetectorElement::g4PhysicalVolume() const {
  return *m_g4physVol;
}

}  // namespace Acts
