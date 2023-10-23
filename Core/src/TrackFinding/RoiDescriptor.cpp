// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// TODO: update to C++17 style
#include "Acts/TrackFinding/RoiDescriptor.hpp"

#include <cmath>
#include <sstream>

namespace Acts {

Acts::RoiDescriptor::RoiDescriptor(double eta, double etaMinus, double etaPlus,
                                   double phi, double phiMinus, double phiPlus,
                                   double zed, double zedMinus, double zedPlus)
    : m_phi(phi),
      m_eta(eta),
      m_zed(zed),
      m_phiMinus(phiMinus),
      m_phiPlus(phiPlus),
      m_etaMinus(etaMinus),
      m_etaPlus(etaPlus),
      m_zedMinus(zedMinus),
      m_zedPlus(zedPlus) {}

Acts::RoiDescriptor::~RoiDescriptor() = default;

}  // namespace Acts
