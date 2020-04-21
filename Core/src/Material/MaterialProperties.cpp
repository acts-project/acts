// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialProperties.hpp"

#include <climits>
#include <limits>
#include <ostream>

static constexpr auto eps = 2 * std::numeric_limits<float>::epsilon();

Acts::MaterialProperties::MaterialProperties(float thickness)
    : m_thickness(thickness) {}

Acts::MaterialProperties::MaterialProperties(float X0, float L0, float Ar,
                                             float Z, float rho,
                                             float thickness)
    : m_material(X0, L0, Ar, Z, rho),
      m_thickness(thickness),
      m_thicknessInX0((X0 > eps) ? (thickness / X0) : 0),
      m_thicknessInL0((L0 > eps) ? (thickness / L0) : 0) {}

Acts::MaterialProperties::MaterialProperties(const Material& material,
                                             float thickness)
    : m_material(material),
      m_thickness(thickness),
      m_thicknessInX0((material.X0() > eps) ? (thickness / material.X0()) : 0),
      m_thicknessInL0((material.L0() > eps) ? (thickness / material.L0()) : 0) {
}

Acts::MaterialProperties::MaterialProperties(
    const std::vector<MaterialProperties>& layers)
    : MaterialProperties() {
  // use double for computations to avoid precision loss
  double Ar = 0.0;
  double Z = 0.0;
  double weight = 0.0;
  double thickness = 0.0;
  double thicknessInX0 = 0.0;
  double thicknessInL0 = 0.0;
  // sum-up contributions from each layer
  for (const auto& layer : layers) {
    const auto& mat = layer.material();
    // weight of the layer assuming a unit area, i.e. volume = thickness*1*1
    const auto layerWeight = mat.massDensity() * layer.thickness();
    // Ar,Z are weighted by mass
    Ar += mat.Ar() * layerWeight;
    Z += mat.Z() * layerWeight;
    weight += layerWeight;
    // thickness and relative thickness in X0,L0 are strictly additive
    thickness += layer.thickness();
    thicknessInX0 += layer.thicknessInX0();
    thicknessInL0 += layer.thicknessInL0();
  }
  // store averaged material constants
  Ar /= weight;
  Z /= weight;
  // this is weight/volume w/ volume = thickness*unitArea = thickness*1*1
  const auto density = weight / thickness;
  const auto X0 = thickness / thicknessInX0;
  const auto L0 = thickness / thicknessInL0;
  m_material = Material(X0, L0, Ar, Z, density);
  // thickness properties do not need to be averaged
  m_thickness = thickness;
  m_thicknessInX0 = thicknessInX0;
  m_thicknessInL0 = thicknessInL0;
}

void Acts::MaterialProperties::scaleThickness(float scale) {
  m_thickness *= scale;
  m_thicknessInX0 *= scale;
  m_thicknessInL0 *= scale;
}

std::ostream& Acts::operator<<(std::ostream& os,
                               const MaterialProperties& materialProperties) {
  os << materialProperties.material()
     << "|t=" << materialProperties.thickness();
  return os;
}
