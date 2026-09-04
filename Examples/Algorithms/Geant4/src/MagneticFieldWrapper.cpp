// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Geant4/UnitConversion.hpp"

#include <utility>

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

namespace ActsExamples::Geant4 {

MagneticFieldWrapper::MagneticFieldWrapper(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4MagneticField(), m_cfg(cfg), m_logger(std::move(logger)) {}

void MagneticFieldWrapper::GetFieldValue(const G4double point[4],
                                         G4double* bField) const {
  assert(bField != nullptr);

  auto bCache = m_cfg.magneticField->makeCache(Acts::MagneticFieldContext());

  auto fieldRes = m_cfg.magneticField->getField(
      {convertLengthToActs * point[0], convertLengthToActs * point[1],
       convertLengthToActs * point[2]},
      bCache);
  if (!fieldRes.ok()) {
    ACTS_ERROR("Field lookup error: " << fieldRes.error());
    return;
  }
  // Get the field now
  const Acts::Vector3& field = *fieldRes;

  bField[0] = convertFieldToGeant4 * field[0];
  bField[1] = convertFieldToGeant4 * field[1];
  bField[2] = convertFieldToGeant4 * field[2];
}

}  // namespace ActsExamples::Geant4
