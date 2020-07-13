// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Material.hpp"

#include <cmath>
#include <ostream>

#include "Acts/Utilities/Units.hpp"

namespace {
enum MaterialClassificationNumberIndices {
  eRadiationLength = 0,
  eInteractionLength = 1,
  eRelativeAtomicMass = 2,
  eNuclearCharge = 3,
  eMassDensity = 4,
};
}

Acts::Material::Material(const ActsVectorF<5>& parameters)
    : Material(parameters[eRadiationLength], parameters[eInteractionLength],
               parameters[eRelativeAtomicMass], parameters[eNuclearCharge],
               parameters[eMassDensity]) {}

float Acts::Material::molarElectronDensity() const {
  using namespace Acts::UnitLiterals;

  // no material, no electron density
  if (!(*this)) {
    return 0.0f;
  }

  // Avogadro constant
  constexpr double Na = 6.02214076e23;
  // perform computations in double precision. due to the native units we might
  // and up with ratios of large numbers and need to keep precision.
  // convert relativ atomic mass Ar to atom mass in native units
  const double atomicMass = m_ar * 1_u;
  // compute molar atomic density
  //   [mass density / atom mass / Avogadro constant]
  //   [mass density / (atom mass * Avogadro constant)]
  // = (mass / volume) / (mass * (1/mol))
  // = (mass * mol) / (mass * volume)
  const double molarAtomicDensity =
      static_cast<double>(m_rho) / atomicMass / Na;
  // each atom has Z electrons
  return static_cast<float>(m_z * molarAtomicDensity);
}

float Acts::Material::meanExcitationEnergy() const {
  using namespace Acts::UnitLiterals;

  // use approximative computation as defined in ATL-SOFT-PUB-2008-003
  return 16_eV * std::pow(m_z, 0.9f);
}

Acts::ActsVectorF<5> Acts::Material::classificationNumbers() const {
  ActsVectorF<5> parameters;
  parameters[eRadiationLength] = m_x0;
  parameters[eInteractionLength] = m_l0;
  parameters[eRelativeAtomicMass] = m_ar;
  parameters[eNuclearCharge] = m_z;
  parameters[eMassDensity] = m_rho;
  return parameters;
}

std::ostream& Acts::operator<<(std::ostream& os, const Material& material) {
  if (!material) {
    os << "vacuum";
  } else {
    os << "X0=" << material.X0();
    os << "|L0=" << material.L0();
    os << "|Ar=" << material.Ar();
    os << "|Z=" << material.Z();
    os << "|rho=" << material.massDensity();
  }
  return os;
}
