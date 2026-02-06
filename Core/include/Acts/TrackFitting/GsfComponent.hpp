// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

namespace Acts {

/// Encapsulates a component of a Gaussian mixture as used by the GSF
struct GsfComponent {
  /// Weight of this component in the Gaussian mixture
  double weight = 0;
  /// Bound track parameters for this component
  BoundVector boundPars = BoundVector::Zero();
  /// Covariance matrix for the bound track parameters
  BoundMatrix boundCov = BoundMatrix::Zero();
};

}  // namespace Acts
