// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/AlignedGeant4DetectorElement.hpp"

std::shared_ptr<ActsExamples::AlignedGeant4DetectorElement>
ActsExamples::alignedGeant4DetectorElementFactory(
    std::shared_ptr<Acts::Surface> surface, const G4VPhysicalVolume& g4physVol,
    const Acts::Transform3& toGlobal, double thickness) {
  return std::make_shared<AlignedGeant4DetectorElement>(
      std::move(surface), g4physVol, toGlobal, thickness);
}
