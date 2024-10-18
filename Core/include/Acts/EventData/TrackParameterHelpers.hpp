// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <cmath>
#include <numbers>

namespace Acts {

/// Check if a bound vector is valid. This checks the following:
/// - All values are finite
/// - The phi value is in the range [-pi, pi)
/// - The theta value is in the range [0, pi]
/// - The q/p value is non-zero
///
/// @param v The bound vector to check
/// @param epsilon The epsilon to use for the checks
///
/// @return True if the bound vector is valid
bool isBoundVectorValid(const BoundVector& v, double epsilon = 1e-6);

/// Check if a free vector is valid. This checks the following:
/// - All values are finite
/// - Direction is normalized
/// - The q/p value is non-zero
///
/// @param v The free vector to check
/// @param epsilon The epsilon to use for the checks
///
/// @return True if the free vector is valid
bool isFreeVectorValid(const FreeVector& v, double epsilon = 1e-6);

}  // namespace Acts
