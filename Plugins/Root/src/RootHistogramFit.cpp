// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include "Acts/Utilities/GaussianHistogramFitError.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <memory>

#include <TFitResult.h>
#include <TH1.h>

using Acts::Experimental::GaussianHistogramFitError;
using Acts::Experimental::GaussianHistogramFitResult;

namespace ActsPlugins {

Acts::Result<GaussianHistogramFitResult> RootHistogramFit::fit(
    const Acts::Experimental::Histogram1& hist) const {
  const std::unique_ptr<TH1F> rootHist = toRoot(hist);
  const TFitResultPtr result = rootHist->Fit("gaus", "SQ0", nullptr);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return Acts::Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::FitFailed);
  }

  return GaussianHistogramFitResult{result->Parameter(1), result->Parameter(2),
                                    result->ParError(1), result->ParError(2)};
}

Acts::Result<GaussianHistogramFitResult> RootHistogramFit::fit(
    const Acts::Experimental::Histogram1& hist, double xMin,
    double xMax) const {
  const std::unique_ptr<TH1F> rootHist = toRoot(hist);
  const TFitResultPtr result =
      rootHist->Fit("gaus", "SRQ0", nullptr, xMin, xMax);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return Acts::Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::FitFailed);
  }

  return GaussianHistogramFitResult{result->Parameter(1), result->Parameter(2),
                                    result->ParError(1), result->ParError(2)};
}

}  // namespace ActsPlugins
