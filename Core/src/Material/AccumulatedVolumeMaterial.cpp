// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"

void Acts::AccumulatedVolumeMaterial::accumulate(const Material& mat) {
  // If nothing is set it is vacuum
  if (!mat) {
    m_vacuumEntries++;
  } else {
    // Replace the vacuum by matter or add matter to matter
    if (m_totalX0 == std::numeric_limits<float>::infinity()) {
      m_totalX0 = mat.X0();
    } else {
      m_totalX0 += mat.X0();
    }
    if (m_totalL0 == std::numeric_limits<float>::infinity()) {
      m_totalL0 = mat.L0();
    } else {
      m_totalL0 += mat.L0();
    }
    m_totalAr += mat.Ar();
    m_totalZ += mat.Z();
    m_totalRho += mat.rho();
    m_materialEntries++;
  }
}

Acts::Material Acts::AccumulatedVolumeMaterial::average() {
  // nothing accumulated, material is vacuum
  if (m_materialEntries == 0) {
    return Material();
  }
  /// The following rescaling is a combination of two steps.
  /// 1) All material entries are averaged.
  /// 2) The numbers are rescaled by a material-to-totalEntries factor. This
  /// rescaling is performed by dividing A, Z and rho by this factor and
  /// multiplying X0 and L0 by it.
  float scalor = m_materialEntries * m_materialEntries;
  float totalEntries = (float)(m_vacuumEntries + m_materialEntries);
  // Create the material
  return Material(m_totalX0 * totalEntries / scalor,
                  m_totalL0 * totalEntries / scalor, m_totalAr / totalEntries,
                  m_totalZ / totalEntries, m_totalRho / totalEntries);
}
