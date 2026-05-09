// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include "CLHEP/Units/SystemOfUnits.h"

namespace ActsExamples::Geant4 {

constexpr double convertLengthToActs = Acts::UnitConstants::mm / CLHEP::mm;
constexpr double convertTimeToActs = Acts::UnitConstants::ns / CLHEP::ns;
constexpr double convertEnergyToActs = Acts::UnitConstants::GeV / CLHEP::GeV;
constexpr double convertFieldToActs = Acts::UnitConstants::T / CLHEP::tesla;
constexpr double convertDensityToActs =
    (Acts::UnitConstants::g / Acts::UnitConstants::mm3) /
    (CLHEP::gram / CLHEP::mm3);

constexpr double convertLengthToGeant4 = CLHEP::mm / Acts::UnitConstants::mm;
constexpr double convertTimeToGeant4 = CLHEP::ns / Acts::UnitConstants::ns;
constexpr double convertEnergyToGeant4 = CLHEP::GeV / Acts::UnitConstants::GeV;
constexpr double convertFieldToGeant4 = CLHEP::tesla / Acts::UnitConstants::T;

}  // namespace ActsExamples::Geant4
