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
Acts::AccumulatedVolumeMaterial::totalAverage()
{
  if (m_totalEntries > 0) {
    double totalScalor = 1. / (double)m_totalEntries;
    m_totalX0 *= totalScalor;
    m_totalL0 *= totalScalor;
    m_totalA *= totalScalor;
    m_totalZ *= totalScalor;
    m_totalRho *= totalScalor;
    // Create the material
    return Material(m_totalX0, m_totalL0, m_totalA, m_totalZ, m_totalRho);
  }
  // Create vacuum
  return Material();
}