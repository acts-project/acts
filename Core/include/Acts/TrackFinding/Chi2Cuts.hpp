// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <array>

namespace Acts {

/// Chi2 cut for each measurement dimension, index 0 is not used
using Chi2CutPerDimension = std::array<double, kMeasurementSizeMax + 1>;

/// Convert a chi2 cut for one-dimensional measurements to the cut with the
/// same p-value for each measurement dimension.
///
/// A cut which is not finite, not positive, or so far in the tail that its
/// p-value is not representable, is the same for every dimension.
///
/// @param cut1D the chi2 cut for one degree of freedom
/// @return the chi2 cut for each number of degrees of freedom
Chi2CutPerDimension chi2CutsWithSameProbability(double cut1D);

}  // namespace Acts
