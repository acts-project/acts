// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialProperties.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/MaterialProperties.hpp"

#include <climits>
#include <ostream>

Acts::MaterialProperties::MaterialProperties(float thickness)
    : m_thickness(thickness) {}

Acts::MaterialProperties::MaterialProperties(float Xo, float Lo, float averageA,
                                             float averageZ, float averageRho,
                                             float thickness)
    : m_material(Xo, Lo, averageA, averageZ, averageRho),
      m_thickness(thickness),
      m_dInX0(Xo * Xo > 10e-10 ? thickness / Xo : 0.),
      m_dInL0(Lo * Lo > 10e-10 ? thickness / Lo : 0.) {}

Acts::MaterialProperties::MaterialProperties(const Material& material,
                                             float thickness)
    : m_material(material),
      m_thickness(thickness),
      m_dInX0(material.X0() * material.X0() > 10e-10 ? thickness / material.X0()
                                                     : 0.),
      m_dInL0(material.L0() * material.L0() > 10e-10 ? thickness / material.L0()
                                                     : 0.) {}

Acts::MaterialProperties::MaterialProperties(
    const std::vector<MaterialProperties>& matLayers, bool normalize)
    : MaterialProperties() {
  double rho = 0.;
  double A = 0.;
  double Z = 0.;
  double X0 = 0.;
  double L0 = 0.;

  for (auto& mat : matLayers) {
    // thickness in X0 and L0 are strictly additive
    m_dInX0 += mat.thicknessInX0();
    m_dInL0 += mat.thicknessInL0();
    double t = mat.thickness();
    double r = mat.material().rho();
    m_thickness += t;
    // density scales with thickness
    rho += r * t;
    // A/Z scale with thickness * density
    A += mat.material().A() * r * t;
    Z += mat.material().Z() * r * t;
  }
  // Now create the average
  X0 = m_thickness / m_dInX0;
  L0 = m_thickness / m_dInL0;
  A /= rho;
  Z /= rho;
  rho /= m_thickness;
  // set the material
  m_material = Material(X0, L0, A, Z, rho);
  if (normalize) {
    scaleToUnitThickness();
  }
}

Acts::MaterialProperties& Acts::MaterialProperties::operator*=(float scale) {
  // assuming rescaling of the material thickness
  m_dInX0 *= scale;
  m_dInL0 *= scale;
  m_thickness *= scale;
  return (*this);
}

void Acts::MaterialProperties::scaleToUnitThickness() {
  // And 'condense to unit thickness' if configured
  double t = thickness();
  double X0 = m_material.X0() / t;
  double L0 = m_material.L0() / t;
  double A = m_material.A();
  double Z = m_material.Z();
  double rho = m_material.rho() * t;
  m_material = Material(X0, L0, A, Z, rho);
  m_thickness = 1.;
}

std::ostream& Acts::operator<<(std::ostream& os,
                               const MaterialProperties& materialProperties) {
  os << materialProperties.material()
     << "|t=" << materialProperties.thickness();
  return os;
}
