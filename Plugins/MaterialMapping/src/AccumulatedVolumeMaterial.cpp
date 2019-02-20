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
  m_eventX0 += mat.X0();
  m_eventL0 += mat.L0();
  m_eventA += mat.A();
  m_eventZ += mat.Z();
  m_eventRho += mat.rho();
  m_totalEntries++;
}

Acts::Material
Acts::AccumulatedVolumeMaterial::eventAverage()
{
  if (m_totalEntires > 0) {
    double eventScalor = 1. / (double)m_totalEntries;
    m_totalX0 *= eventScalor;
    m_totalL0 *= eventScalor;
    m_totalA *= eventScalor;
    m_totalZ *= eventScalor;
    m_totalRho *= eventScalor;
    // Create the material
    return Material(m_totalX0, m_totalL0, m_totalA, m_totalZ, m_totalRho);
  } else {
    // Create vacuum
    return Material();
  }
}