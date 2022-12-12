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

#include "G4VPhysicalVolume.hh"

namespace Acts {

class IGeant4PhysicalVolumeSelector {
 public:
  virtual ~IGeant4PhysicalVolumeSelector(){};
};

/// The selector delegate for G4VPhysical volumes
using Geant4PhysicalVolumeSelector =
    OwningDelegate<bool(const G4VPhysicalVolume&),
                   IGeant4PhysicalVolumeSelector>;

namespace Geant4PhysicalVolumeSelectors {

struct AllSelector : public IGeant4PhysicalVolumeSelector {
  bool select(const G4VPhysicalVolume& /*unused*/) const { return true; }
};

static inline Geant4PhysicalVolumeSelector generateAllSelector() {
  Geant4PhysicalVolumeSelector allSelector;
  auto all = std::make_unique<const AllSelector>();
  allSelector.connect<&AllSelector::select>(std::move(all));
  return allSelector;
}

struct NameSelector : public IGeant4PhysicalVolumeSelector {
  std::string name = "";
  bool exact = false;

  NameSelector(const std::string& n, bool e = false) : name(n), exact(e) {}

  bool select(const G4VPhysicalVolume& g4PhysVol) const {
    std::string volumeName = g4PhysVol.GetName();
    return exact ? (volumeName == name)
                 : volumeName.find(name) != std::string::npos;
  }
};

static inline Geant4PhysicalVolumeSelector generateNameSelector(
    const std::string& name, bool exact = false) {
  Geant4PhysicalVolumeSelector nameSelector;
  auto nsel = std::make_unique<const NameSelector>(name, exact);
  nameSelector.connect<&NameSelector::select>(std::move(nsel));
  return nameSelector;
}

/*
struct ExtentSelector : public IGeant4PhysicalVolumeSelector {
    Extent restriction;

  bool select(const G4VPhysicalVolume& g4PhysVol) const {

  }

};
*/

}  // namespace Geant4PhysicalVolumeSelectors
}  // namespace Acts
