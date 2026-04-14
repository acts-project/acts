// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Geant4/Geant4DetectorElement.hpp"

#include "Acts/Surfaces/Surface.hpp"


using namespace Acts;

namespace ActsPlugins {

Geant4DetectorElement::Geant4DetectorElement(const G4VPhysicalVolume& g4physVol,
                                             const Transform3& toGlobal,
                                             double thickness)
    : m_g4physVol(&g4physVol), m_toGlobal(toGlobal), m_thickness(thickness) {}

const Transform3& Geant4DetectorElement::localToGlobalTransform(
    const GeometryContext& /*gctx*/) const {
  return m_toGlobal;
}

const Surface& Geant4DetectorElement::surface() const {
  return *SurfacePlacementBase::surface();
}

Surface& Geant4DetectorElement::surface() {
  return const_cast<Surface&>(*SurfacePlacementBase::surface());
}

double Geant4DetectorElement::thickness() const {
  return m_thickness;
}

const G4VPhysicalVolume& Geant4DetectorElement::g4PhysicalVolume() const {
  return *m_g4physVol;
}

}  // namespace ActsPlugins
