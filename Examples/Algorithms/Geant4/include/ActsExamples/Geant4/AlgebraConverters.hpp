// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "G4ThreeVector.hh"

namespace ActsExamples::Geant4 {
/// @brief Converts a G4 3-vector position into an Acts 3-vector
Acts::Vector3 convertPosition(const G4ThreeVector& g4vec);
/// @brief Converts a G4 4-vector position into an Acts 4-vector
Acts::Vector4 convertPosition(const G4ThreeVector& g4vec, const double time);
/// @brief Converts a G4 momentum vector into an Acts momentum vector
Acts::Vector4 convertMomentum(const G4ThreeVector& g4vec, const double energy);
/// @brief Converts a G4 direction vector into an Acts direction vector
Acts::Vector3 convertDirection(const G4ThreeVector& g4vec);
/// @brief Converts a Acts 3-vector position into a G4 3-vector
G4ThreeVector convertPosition(const Acts::Vector3& actsVec);
/// @brief Converts a Acts 3-direction position into a G4 3-direction
G4ThreeVector convertDirection(const Acts::Vector3& actsVec);

}  // namespace ActsExamples::Geant4
