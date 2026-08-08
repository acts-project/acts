// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/NumericalMinimizationError.hpp"

#include <string>

namespace {

class NumericalMinimizationErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final {
    return "NumericalMinimizationError";
  }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::detail::NumericalMinimizationError;

    switch (static_cast<NumericalMinimizationError>(c)) {
      case NumericalMinimizationError::DidNotConverge:
        return "Nelder-Mead simplex search did not converge";
      case NumericalMinimizationError::NoFiniteVertex:
        return "No vertex of the initial simplex evaluated to a finite value";
      case NumericalMinimizationError::InvalidStep:
        return "A finite-difference step was not strictly positive";
      case NumericalMinimizationError::NotPositiveDefinite:
        return "The finite-difference Hessian was not positive definite";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::detail::make_error_code(
    Acts::detail::NumericalMinimizationError e) {
  static NumericalMinimizationErrorCategory c;
  return {static_cast<int>(e), c};
}
