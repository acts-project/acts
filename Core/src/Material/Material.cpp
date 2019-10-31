// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Material.hpp"

#include <cmath>
#include <ostream>

#include "Acts/Utilities/Units.hpp"

Acts::Material::Material(float iX0, float iL0, float iA, float iZ, float iRho)
    : m_x0(iX0), m_l0(iL0), m_a(iA), m_z(iZ), m_rho(iRho) {}

Acts::Material::Material(const ActsVectorF<5>& parameters)
    : Material(parameters[eX0], parameters[eL0], parameters[eA], parameters[eZ],
               parameters[eRho]) {}

float Acts::Material::electronDensity() const {
  // no material, no electron density
  if (!(*this)) {
    return 0.0f;
  }

  // atomic mass constant to convert from relative atomic mass to native unit.
  constexpr double Da = 931.49410242_MeV;
  // Avogadro constant
  constexpr double Na = 6.02214076e23;

  // perform computations in double precision. due to our choice of native
  // units we might and up with a ratio of large numbers and need to
  // avoid precision loss.
  const double A = m_a;
  const double Z = m_z;
  const double rho = m_rho;
  // atomic mass in native units is A * Da
  const auto atomMass = (A * Da);
  // each atom has Z electrons, thus electrons per mass is Z / (A * Da)
  // units are  [Z / (A * Da)] = 1 / (1 * mass)
  const auto electronsPerMass = Z / atomMass;
  // electrons are thus given by (Z / (A * Da)) * rho
  // units are [(Z / (A * Da)) * rho] = (1/mass) * mass/volume = 1/volume
  const auto electronDensity = electronsPerMass * rho;
  // convert from numbers of electrons to mol.
  const auto molDensity = electronDensity / Na;

  return molDensity;
}

float Acts::Material::meanExcitationEnergy() const {
  using namespace Acts::UnitLiterals;
  // use approximative computation as defined in ATL-SOFT-PUB-2008-003
  return 16_eV * std::pow(m_z, 0.9f);
}

Acts::ActsVectorF<5> Acts::Material::classificationNumbers() const {
  ActsVectorF<5> parameters;
  parameters[eX0] = m_x0;
  parameters[eL0] = m_l0;
  parameters[eA] = m_a;
  parameters[eZ] = m_z;
  parameters[eRho] = m_rho;
  return parameters;
}

std::ostream& Acts::operator<<(std::ostream& os, const Material& material) {
  if (!material) {
    os << "vacuum";
  } else {
    os << "X0=" << material.X0();
    os << "|L0=" << material.L0();
    os << "|A=" << material.A();
    os << "|Z=" << material.Z();
    os << "|rho=" << material.rho();
  }
  return os;
}
