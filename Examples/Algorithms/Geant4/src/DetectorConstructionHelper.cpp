// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/DetectorConstructionHelper.hpp"

#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeG4DetectorConstruction.hpp"
#ifdef ACTS_GEANT4_DD4HEP
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DDG4/DDG4DetectorConstruction.hpp"
#endif

#include <stdexcept>

G4VUserDetectorConstruction *ActsExamples::getG4DetectorContruction(
    const ActsExamples::IBaseDetector &detector) {
  if (auto telescope =
          dynamic_cast<const ActsExamples::Telescope::TelescopeDetector *>(
              &detector)) {
    return new ActsExamples::Telescope::TelescopeG4DetectorConstruction(
        telescope->config);
#ifdef ACTS_GEANT4_DD4HEP
  } else if (auto dd4hep =
                 dynamic_cast<const ActsExamples::DD4hep::DD4hepDetector *>(
                     &detector)) {
    return new ActsExamples::DDG4DetectorConstruction(
        *dd4hep->geometryService->lcdd());
#endif
  }

  throw std::invalid_argument("given detector has no Geant4 construction");
}
