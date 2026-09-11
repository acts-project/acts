// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsRoiDescriptor.hpp"

#include <cmath>

namespace Acts::Experimental {

GbtsRoiDescriptor::GbtsRoiDescriptor(double etaMin, double etaMax, double zMin,
                                     double zMax)
    : m_zMin(zMin), m_zMax(zMax) {
  // cot(2 * atan(exp(-eta))) is sinh(eta)
  m_dzdrMin = std::sinh(etaMin);
  m_dzdrMax = std::sinh(etaMax);
}

}  // namespace Acts::Experimental
