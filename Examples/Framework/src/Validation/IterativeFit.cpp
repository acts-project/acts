// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/IterativeFit.hpp"

#include <utility>

namespace ActsExamples {

std::optional<HistogramFitResult> iterativeFit(
    const HistogramFitFunction& fitFn,
    const Acts::Experimental::Histogram1& hist, double sigmaRange,
    int iterations, const Acts::Logger& logger) {
  std::optional<HistogramFitResult> result = fitFn(hist, std::nullopt);
  if (!result.has_value()) {
    ACTS_DEBUG("Failed to fit initial Gaussian to '" << hist.name() << "'");
    return result;
  }

  for (int i = 0; i < iterations - 1; ++i) {
    const double mean = std::get<0>(*result);
    const double sigma = std::get<1>(*result);
    const double xMin = mean - sigmaRange * sigma;
    const double xMax = mean + sigmaRange * sigma;

    std::optional<HistogramFitResult> restricted =
        fitFn(hist, HistogramFitRange{xMin, xMax});
    if (!restricted.has_value()) {
      ACTS_DEBUG("Failed to fit iteration " << i << " Gaussian to '"
                                            << hist.name() << "'");
      return restricted;
    }

    result = std::move(restricted);
  }

  return result;
}

}  // namespace ActsExamples
