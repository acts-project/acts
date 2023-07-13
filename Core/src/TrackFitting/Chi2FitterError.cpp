// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/Chi2FitterError.hpp"

#include <string>

namespace {

class Chi2FitterErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "Chi2FitterError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::Experimental::Chi2FitterError;

    switch (static_cast<Chi2FitterError>(c)) {
      case Chi2FitterError::NoMeasurementFound:
        return "No measurement detected during the propagation";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::Experimental::make_error_code(
    Acts::Experimental::Chi2FitterError e) {
  static Chi2FitterErrorCategory c;
  return {static_cast<int>(e), c};
}
