// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/TrackParamsEstimationError.hpp"

#include <string>

namespace {

class TrackParamsEstimationErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final {
    return "TrackParamsEstimationError";
  }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::TrackParamsEstimationError;

    switch (static_cast<TrackParamsEstimationError>(c)) {
      case TrackParamsEstimationError::NotEnoughSpacePoints:
        return "At least three space points are required";
      case TrackParamsEstimationError::DegenerateFit:
        return "The space point configuration leads to a degenerate fit";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::TrackParamsEstimationError e) {
  static TrackParamsEstimationErrorCategory c;
  return {static_cast<int>(e), c};
}
