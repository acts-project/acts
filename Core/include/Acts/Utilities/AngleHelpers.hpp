// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Acts::AngleHelpers {

/// Calculate the pseudorapidity from the polar angle theta.
///
/// @param theta is the polar angle in radian towards the z-axis.
///
/// @return the pseudorapidity towards the z-axis.
template <typename Scalar>
Scalar etaFromTheta(Scalar theta) {
  return -std::log(std::tan(0.5 * theta));
}

/// Calculate the polar angle theta from the pseudorapidity.
///
/// @param eta is the pseudorapidity towards the z-axis.
///
/// @return the polar angle in radian towards the z-axis.
template <typename Scalar>
Scalar thetaFromEta(Scalar eta) {
  return 2 * std::atan(std::exp(-eta));
}

}  // namespace Acts::AngleHelpers
