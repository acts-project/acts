// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Material.hpp"

#include "Acts/Definitions/Units.hpp"

#include <cmath>
#include <ostream>

namespace {
enum MaterialClassificationNumberIndices {
  eRadiationLength = 0,
  eInteractionLength = 1,
  eRelativeAtomicMass = 2,
  eNuclearCharge = 3,
  eMolarDensity = 4,
  eMeanExcitationEnergy = 5,
};
}  // namespace

Acts::Material Acts::Material::fromMassDensity(float x0, float l0, float ar,
                                               float z, float massRho,
                                               std::optional<float> I) {
  using namespace Acts::UnitLiterals;

  Material mat;
  mat.m_x0 = x0;
  mat.m_l0 = l0;
  mat.m_ar = ar;
  mat.m_z = z;
  // mass density is defined as
  //
  //     mass-density = atomic-mass * number-of-atoms / volume
  //                  = atomic-mass * molar-density * avogadro-constant
  // -> molar-density = mass-density / (atomic-mass * avogadro-constant)
  //
  // with the atomic mass given by
  //
  //      atomic-mass = relative-atomic-mass * atomic-mass-unit
  //
  // perform computations in double precision to avoid loss of precision
  const double atomicMass = static_cast<double>(ar) * 1_u;
  mat.m_molarRho =
      static_cast<double>(massRho) / (atomicMass * UnitConstants::Avogadro);

  if (I) {
    mat.m_I = *I;
  } else {
    // use approximative computation as defined in ATL-SOFT-PUB-2008-003
    mat.m_I = 16_eV * std::pow(z, 0.9f);
  }

  return mat;
}

Acts::Material Acts::Material::fromMolarDensity(float x0, float l0, float ar,
                                                float z, float molarRho,
                                                std::optional<float> I) {
  using namespace Acts::UnitLiterals;

  Material mat;
  mat.m_x0 = x0;
  mat.m_l0 = l0;
  mat.m_ar = ar;
  mat.m_z = z;
  mat.m_molarRho = molarRho;

  if (I) {
    mat.m_I = *I;
  } else {
    // use approximative computation as defined in ATL-SOFT-PUB-2008-003
    mat.m_I = 16_eV * std::pow(z, 0.9f);
  }

  return mat;
}

Acts::Material::Material(const ParametersVector& parameters)
    : m_x0(parameters[eRadiationLength]),
      m_l0(parameters[eInteractionLength]),
      m_ar(parameters[eRelativeAtomicMass]),
      m_z(parameters[eNuclearCharge]),
      m_molarRho(parameters[eMolarDensity]) {}

float Acts::Material::massDensity() const {
  using namespace Acts::UnitLiterals;

  // perform computations in double precision to avoid loss of precision
  const double atomicMass = static_cast<double>(m_ar) * 1_u;
  const double numberDensity =
      static_cast<double>(m_molarRho) * UnitConstants::Avogadro;
  return atomicMass * numberDensity;
}

Acts::Material::ParametersVector Acts::Material::parameters() const {
  ParametersVector parameters;
  parameters[eRadiationLength] = m_x0;
  parameters[eInteractionLength] = m_l0;
  parameters[eRelativeAtomicMass] = m_ar;
  parameters[eNuclearCharge] = m_z;
  parameters[eMolarDensity] = m_molarRho;
  parameters[eMeanExcitationEnergy] = m_I;
  return parameters;
}

std::ostream& Acts::operator<<(std::ostream& os, const Material& material) {
  if (!material) {
    os << "vacuum";
  } else {
    os << "x0=" << material.X0();
    os << "|l0=" << material.L0();
    os << "|ar=" << material.Ar();
    os << "|z=" << material.Z();
    os << "|molar_rho=" << material.molarDensity();
    os << "|I=" << material.meanExcitationEnergy();
  }
  return os;
}
