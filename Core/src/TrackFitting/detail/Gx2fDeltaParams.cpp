// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// GX2F backend. One heavy Eigen decomposition per translation unit to keep
// the peak compile memory low.

#include "Acts/TrackFitting/GlobalChiSquareFitterBackend.hpp"

Eigen::VectorXd Acts::Experimental::computeGx2fDeltaParams(
    const Acts::Experimental::Gx2fSystem& extendedSystem) {
  return extendedSystem.aMatrix().colPivHouseholderQr().solve(
      extendedSystem.bVector());
}
