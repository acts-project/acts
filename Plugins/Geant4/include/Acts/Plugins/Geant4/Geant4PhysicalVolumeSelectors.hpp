// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Delegate.hpp"

#include <string>
#include <vector>

#include "G4VPhysicalVolume.hh"

namespace Acts {

/// Interface class for selectors from physical volumes
class IGeant4PhysicalVolumeSelector {
 public:
  virtual ~IGeant4PhysicalVolumeSelector() = default;
  /// @brief  The main interface method
  /// @param g4Phys the physical volume to be checked
  /// @return a boolen indicating if it should be selected or not
  virtual bool select(const G4VPhysicalVolume& g4Phys) const = 0;
};

namespace Geant4PhysicalVolumeSelectors {

/// @brief  Struct that selects all G4VPhysicalVolume objects
struct AllSelector : public IGeant4PhysicalVolumeSelector {
  bool select(const G4VPhysicalVolume& /*unused*/) const final { return true; }
};

/// @brief Struct that selects G4VPhysicalVolume objects
/// that match one of the provided names, exact or partially
struct NameSelector : public IGeant4PhysicalVolumeSelector {
  std::vector<std::string> names = {};
  bool exact = false;

  /// Constructor with arguments
  /// @param ns the provided list of names
  /// @param e whether to select them exact or not
  NameSelector(const std::vector<std::string>& ns, bool e = false)
      : names(ns), exact(e) {}

  /// Secect function for the volume
  /// @param g4PhysVol the volume that is checked
  /// @return a boolean indicating the selection
  bool select(const G4VPhysicalVolume& g4PhysVol) const final {
    std::string volumeName = g4PhysVol.GetName();
    bool matched = false;
    for (const auto& name : names) {
      matched = exact ? (volumeName == name)
                      : volumeName.find(name) != std::string::npos;
      if (matched) {
        break;
      }
    }
    return matched;
  }
};

}  // namespace Geant4PhysicalVolumeSelectors
}  // namespace Acts
