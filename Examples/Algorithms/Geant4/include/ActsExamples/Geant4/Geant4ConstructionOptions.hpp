// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>

namespace ActsExamples::Geant4 {
class RegionCreator;
}  // namespace ActsExamples::Geant4

namespace ActsExamples {

/// @brief Helper struct to pass additional options for the G4 Geometry construction.
struct Geant4ConstructionOptions {
  std::vector<std::shared_ptr<Geant4::RegionCreator>> regionCreators;
};

}  // namespace ActsExamples
