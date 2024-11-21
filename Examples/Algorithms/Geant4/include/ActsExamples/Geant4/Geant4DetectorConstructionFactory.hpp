// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <G4VUserDetectorConstruction.hh>

namespace ActsExamples {

namespace Geant4 {
class RegionCreator;
}

class Geant4DetectorConstructionFactory {
 public:
  virtual ~Geant4DetectorConstructionFactory() = default;

  virtual std::unique_ptr<G4VUserDetectorConstruction> factorize(
      const std::vector<std::shared_ptr<Geant4::RegionCreator>>& regionCreators)
      const = 0;
};

}  // namespace ActsExamples
