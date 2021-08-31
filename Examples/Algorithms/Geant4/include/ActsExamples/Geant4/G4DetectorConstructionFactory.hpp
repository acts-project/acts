// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "G4VUserDetectorConstruction.hh"

namespace ActsExamples {

/// @brief Class to construct a G4DetectorConstruction
///
class G4DetectorConstructionFactory {
 public:
  /// Virtual destructor for memory safety.
  virtual ~G4DetectorConstructionFactory() = default;

  /// The main virtual method to create a detector construction
  /// @return Unique pointer to the detector construction
  virtual std::unique_ptr<G4VUserDetectorConstruction> operator()() const = 0;
};
}  // namespace ActsExamples