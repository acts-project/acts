// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <memory>

#include <TFitResult.h>
#include <TH1.h>

namespace ActsPlugins {

std::optional<RootHistogramFit::Result> RootHistogramFit::operator()(
    const Acts::Experimental::Histogram1& hist,
    std::optional<Range> range) const {
  const std::unique_ptr<TH1F> rootHist = toRoot(hist);

  const TFitResultPtr result =
      range.has_value()
          ? rootHist->Fit("gaus", (m_config.fitOptions + "R").c_str(), nullptr,
                          range->first, range->second)
          : rootHist->Fit("gaus", m_config.fitOptions.c_str(), nullptr);

  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return std::nullopt;
  }

  return Result{result->Parameter(1), result->Parameter(2), result->ParError(1),
                result->ParError(2)};
}

}  // namespace ActsPlugins
