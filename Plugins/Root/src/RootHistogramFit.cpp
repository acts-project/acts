// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include "Acts/Utilities/GaussianFitError.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <memory>

#include <TFitResult.h>
#include <TH1.h>

using Acts::Experimental::GaussianFitError;
using Acts::Experimental::GaussianFitResult;

namespace ActsPlugins {

Acts::Result<GaussianFitResult> RootHistogramFit::fit(
    const Acts::Experimental::Histogram1& hist) const {
  const std::unique_ptr<TH1F> rootHist = toRoot(hist);
  const TFitResultPtr result = rootHist->Fit("gaus", "SQ0", nullptr);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return Acts::Result<GaussianFitResult>::failure(
        GaussianFitError::FitFailed);
  }

  return GaussianFitResult{result->Parameter(1), result->Parameter(2),
                           result->ParError(1), result->ParError(2)};
}

Acts::Result<GaussianFitResult> RootHistogramFit::fit(
    const Acts::Experimental::Histogram1& hist, double xMin,
    double xMax) const {
  const std::unique_ptr<TH1F> rootHist = toRoot(hist);
  const TFitResultPtr result =
      rootHist->Fit("gaus", "SRQ0", nullptr, xMin, xMax);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return Acts::Result<GaussianFitResult>::failure(
        GaussianFitError::FitFailed);
  }

  return GaussianFitResult{result->Parameter(1), result->Parameter(2),
                           result->ParError(1), result->ParError(2)};
}

}  // namespace ActsPlugins
