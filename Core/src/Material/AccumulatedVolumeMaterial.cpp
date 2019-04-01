// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedVolumeMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"

void
Acts::AccumulatedVolumeMaterial::accumulate(const Material& mat)
{
  // If nothing is set it is vacuum
  if (mat.A() == 0. || mat.Z() == 0. || mat.rho() == 0.) {
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
    m_totalA += mat.A();
    m_totalZ += mat.Z();
    m_totalRho += mat.rho();
    m_materialEntries++;
  }
}

Acts::Material
Acts::AccumulatedVolumeMaterial::average()
{
  if (m_materialEntries > 0) {
    /// The following rescaling is a combination of two steps.
    /// 1) All material entries are averaged.
    /// 2) The numbers are rescaled by a material-to-totalEntries factor. This
    /// rescaling is performed by dividing A, Z and rho by this factor and
    /// multiplying X0 and L0 by it.
    float scalor       = m_materialEntries * m_materialEntries;
    float totalEntries = (float)(m_vacuumEntries + m_materialEntries);
    // Create the material
    return Material(m_totalX0 * totalEntries / scalor,
                    m_totalL0 * totalEntries / scalor,
                    m_totalA / totalEntries,
                    m_totalZ / totalEntries,
                    m_totalRho / totalEntries);
  }
  // Create vacuum
  return Material();
}