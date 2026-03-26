// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GbtsRoiDescriptor.hpp"

#include <cmath>

namespace Acts::Experimental {

GbtsRoiDescriptor::GbtsRoiDescriptor(double eta, double etaMin, double etaMax,
                                     double phi, double phiMin, double phiMax,
                                     double z, double zMin, double zMax)
    : m_phi(phi),
      m_eta(eta),
      m_z(z),
      m_phiMin(phiMin),
      m_phiMax(phiMax),
      m_etaMin(etaMin),
      m_etaMax(etaMax),
      m_zMin(zMin),
      m_zMax(zMax) {
  m_drdzMin = std::tan(2 * std::atan(std::exp(-m_etaMin)));
  m_drdzMax = std::tan(2 * std::atan(std::exp(-m_etaMax)));

  m_dzdrMin = 1 / m_drdzMin;
  m_dzdrMax = 1 / m_drdzMax;
}

}  // namespace Acts::Experimental
