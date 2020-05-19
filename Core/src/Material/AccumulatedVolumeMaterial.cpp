// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include <iostream>

void Acts::AccumulatedVolumeMaterial::accumulate(
    const MaterialProperties& mat) {
  // Replace the vacuum by matter or add matter to matter
  if (m_totalX0 == std::numeric_limits<float>::infinity()) {
    m_totalX0 = mat.thickness() / mat.material().X0();
  } else {
    m_totalX0 += mat.thickness() / mat.material().X0();
  }
  if (m_totalL0 == std::numeric_limits<float>::infinity()) {
    m_totalL0 = mat.thickness() / mat.material().L0();
  } else {
    m_totalL0 += mat.thickness() / mat.material().L0();
  }
  m_totalAr += mat.material().Ar();
  m_totalZ += mat.material().Z();
  m_totalRho += mat.material().massDensity();
  if (m_materialEntries == 0) {
    m_thickness = mat.thickness();
  } else {
    m_thickness += mat.thickness();
  }
  m_materialEntries++;
}

Acts::Material Acts::AccumulatedVolumeMaterial::average() {
  // nothing accumulated, material is vacuum
  if (m_materialEntries == 0) {
    return Material();
  }

  // Create the material
  return Material(m_thickness / m_totalX0, m_thickness / m_totalL0,
                  m_totalAr / m_materialEntries, m_totalZ / m_materialEntries,
                  m_totalRho / m_materialEntries);
}
