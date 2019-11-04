// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

namespace {
/// Scale material properties to unit thickness.
///
/// A helper method to allows to scale a material property for
/// unphysical/blended material to a unit thickness of 1. This is safe for
/// energy loss and multiple scattering application in the material integration
///
/// Scaling to unit thickness changes only X0, L0, and rho.
Acts::Material makeUnitThicknessMaterial(double X0, double L0, double A,
                                         double Z, double density,
                                         double thickness) {
  return {static_cast<float>(X0 / thickness),
          static_cast<float>(L0 / thickness), static_cast<float>(A),
          static_cast<float>(Z), static_cast<float>(density * thickness)};
}
}  // namespace

Acts::MaterialProperties::MaterialProperties(
    const std::vector<MaterialProperties>& layers, bool normalize)
    : MaterialProperties() {
  // use double for computations to avoid precision loss
  double A = 0.;
  double Z = 0.;
  double density = 0.;
  double thickness = 0.;
  for (const auto& layer : layers) {
    const auto& mat = layer.material();
    // A/Z scale with thickness * density
    A += mat.A() * mat.rho() * layer.thickness();
    Z += mat.Z() * mat.rho() * layer.thickness();
    // density scales with thickness
    density += mat.rho() * layer.thickness();
    thickness += layer.thickness();
    // relative thickness in X0 and L0 are strictly additive
    m_dInX0 += layer.thicknessInX0();
    m_dInL0 += layer.thicknessInL0();
  }
  // create the average
  const double X0 = thickness / m_dInX0;
  const double L0 = thickness / m_dInL0;
  A /= density;
  Z /= density;
  density /= thickness;
  // set the material
  if (normalize) {
    m_material = makeUnitThicknessMaterial(X0, L0, A, Z, density, thickness);
    m_thickness = 1.0f;
  } else {
    m_material = Material(X0, L0, A, Z, density);
    m_thickness = thickness;
  }
}

Acts::MaterialProperties& Acts::MaterialProperties::operator*=(float scale) {
  // assuming rescaling of the material thickness
  m_dInX0 *= scale;
  m_dInL0 *= scale;
  m_thickness *= scale;
  return (*this);
}

std::ostream& Acts::operator<<(std::ostream& os,
                               const MaterialProperties& materialProperties) {
  os << materialProperties.material()
     << "|t=" << materialProperties.thickness();
  return os;
}
