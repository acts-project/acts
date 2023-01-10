// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"

Acts::Geant4DetectorElement::Geant4DetectorElement(
    std::shared_ptr<Acts::Surface> surface, Acts::ActsScalar thickness,
    const G4VPhysicalVolume& g4physVol)
    : m_surface(std::move(surface)),
      m_thickness(thickness),
      m_g4physVol(&g4physVol) {}

const Acts::Transform3& Acts::Geant4DetectorElement::transform(
    const GeometryContext& gctx) const {
  return m_surface->transform(gctx);
}

const Acts::Surface& Acts::Geant4DetectorElement::surface() const {
  return (*m_surface);
}

Acts::ActsScalar Acts::Geant4DetectorElement::thickness() const {
  return m_thickness;
}

const G4VPhysicalVolume& Acts::Geant4DetectorElement::g4PhysicalVolume() const {
  return (*m_g4physVol);
}
