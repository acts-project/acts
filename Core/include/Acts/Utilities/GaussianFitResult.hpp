// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::Experimental {

/// @brief Outcome of a Gaussian fit to a 1D histogram
///
/// @note Kept in its own header, separate from @c GaussianFit.hpp, so that
///       generic fit-loop helpers (see @c Acts::detail::iterativeGaussianFit)
///       can depend on the result type without depending on a particular fit
///       implementation.
struct GaussianFitResult {
  /// Fitted mean of the Gaussian
  double mean{};
  /// Fitted standard deviation of the Gaussian
  double sigma{};
  /// Estimated uncertainty on @c mean
  double meanError{};
  /// Estimated uncertainty on @c sigma
  double sigmaError{};
};

}  // namespace Acts::Experimental
