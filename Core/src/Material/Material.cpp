// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/Material.hpp"

#include <ostream>

Acts::Material::Material(float iX0, float iL0, float iA, float iZ, float iRho)
    : m_x0(iX0),
      m_l0(iL0),
      m_a(iA),
      m_z(iZ),
      m_rho(iRho),
      m_ne((0 < iA) ? (iZ / iA) * iRho : 0.0f) {}

Acts::Material::Material(const ActsVectorF<5>& parameters)
    : Material(parameters[eX0], parameters[eL0], parameters[eA], parameters[eZ],
               parameters[eRho]) {}

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
