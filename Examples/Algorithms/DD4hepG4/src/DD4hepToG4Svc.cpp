// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"

#include "ACTFW/Plugins/DD4hepG4/GeoConstruction.hpp"

FW::DD4hepG4::DD4hepToG4Svc::DD4hepToG4Svc(
    const FW::DD4hepG4::DD4hepToG4Svc::Config& cfg)
    : m_cfg(cfg), m_geant4Geometry(nullptr) {}

FW::DD4hepG4::DD4hepToG4Svc::~DD4hepToG4Svc() {
  /// @todo Clarify why this code is commented out
  // delete m_geant4Geometry;
}

FW::ProcessCode FW::DD4hepG4::DD4hepToG4Svc::buildGeant4Geometry() {
  if (m_cfg.dd4hepService->lcdd()) {
    G4VUserDetectorConstruction* detector(
        new FW::DD4hepG4::GeoConstruction(*(m_cfg.dd4hepService->lcdd())));
    m_geant4Geometry = detector;
  }
  /// @todo Clarify why this code is commented out
  // if (!m_geant4Geometry) FW::ProcessCode::ERROR;
  return FW::ProcessCode::SUCCESS;
}

G4VUserDetectorConstruction* FW::DD4hepG4::DD4hepToG4Svc::geant4Geometry() {
  if (!m_geant4Geometry)
    buildGeant4Geometry();
  return m_geant4Geometry;
}
