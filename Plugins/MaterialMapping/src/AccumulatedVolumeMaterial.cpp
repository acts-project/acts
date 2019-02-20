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

#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"

void
Acts::AccumulatedVolumeMaterial::accumulate(const Material& mat)
{
  m_totalX0 += mat.X0();
  m_totalL0 += mat.L0();
  m_totalA += mat.A();
  m_totalZ += mat.Z();
  m_totalRho += mat.rho();
  m_totalEntries++;
}

Acts::Material
Acts::AccumulatedVolumeMaterial::average()
{
  if (m_totalEntries > 0) {
    double scalor = 1. / (double)m_totalEntries;
    // Create the material
    return Material(m_totalX0 * scalor,
                    m_totalL0 * scalor,
                    m_totalA * scalor,
                    m_totalZ * scalor,
                    m_totalRho * scalor);
  }
  // Create vacuum
  return Material();
}