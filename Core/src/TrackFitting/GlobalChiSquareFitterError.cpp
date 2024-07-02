// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitterError.hpp"

#include <string>

namespace {

class GlobalChiSquareFitterErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final {
    return "GlobalChiSquareFitterError";
  }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::Experimental::GlobalChiSquareFitterError;

    switch (static_cast<GlobalChiSquareFitterError>(c)) {
      case GlobalChiSquareFitterError::AIsNotInvertible:
        return "Gx2f: aMatrix is not invertible.";
      case GlobalChiSquareFitterError::DidNotConverge:
        return "Gx2f: Did not converge in 'nUpdateMax' updates.";
      case GlobalChiSquareFitterError::NotEnoughMeasurements:
        return "Gx2f: Not enough measurements.";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::Experimental::make_error_code(
    Acts::Experimental::GlobalChiSquareFitterError e) {
  static GlobalChiSquareFitterErrorCategory c;
  return {static_cast<int>(e), c};
}
