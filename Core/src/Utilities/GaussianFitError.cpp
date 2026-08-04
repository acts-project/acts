// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/GaussianFitError.hpp"

#include <string>

namespace {

class GaussianFitErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "GaussianFitError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::Experimental::GaussianFitError;

    switch (static_cast<GaussianFitError>(c)) {
      case GaussianFitError::EmptyRange:
        return "No bins in the fit range, or their contents summed to zero "
               "or less";
      case GaussianFitError::TooFewNonEmptyBins:
        return "Fewer populated bins than the fitter's configured minimum";
      case GaussianFitError::NoValidSeed:
        return "The initial mean/sigma guess could not be formed";
      case GaussianFitError::NonFiniteParameters:
        return "The minimiser returned a non-finite mean or a non-positive "
               "sigma";
      case GaussianFitError::FitFailed:
        return "The underlying fitter reported failure";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::Experimental::make_error_code(
    Acts::Experimental::GaussianFitError e) {
  static GaussianFitErrorCategory c;
  return {static_cast<int>(e), c};
}
