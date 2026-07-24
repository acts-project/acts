// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"

namespace ActsExamples {

/// Configuration for extracting mean/width profiles from residual/pull
/// histograms via iterative Gaussian fits.
struct ResPlotRefinementConfig {
  /// Minimum number of entries in a bin for it to be included in the
  /// mean/width fit.
  int fitMinEntries = 10;
  /// The range in sigma for the iterative Gaussian fit
  double fitSigmaRange = 3.0;
  /// The maximum number of iterations for the iterative Gaussian fit
  int fitIterations = 3;
  /// Threshold for warning about fit failure fraction in profile extraction.
  double warningThresholdFitFailureFraction = 0.55;
};

/// Write all residual and pull histograms of the given plot tool to the
/// current ROOT directory, together with mean/width profiles extracted via
/// iterative Gaussian fits.
///
/// @param resPlotTool the plot tool holding the filled histograms
/// @param config the profile extraction configuration
/// @param logger logger for fit-failure warnings
void writeResPlots(const ResPlotTool& resPlotTool,
                   const ResPlotRefinementConfig& config,
                   const Acts::Logger& logger);

}  // namespace ActsExamples
