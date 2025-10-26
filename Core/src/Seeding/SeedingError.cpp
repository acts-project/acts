// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedingError.hpp"

#include <string>

namespace {

class SeedingErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "SeedingError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::SeedingError;

    switch (static_cast<SeedingError>(c)) {
      case Acts::SeedingError::InvalidSpacePointsForEstimate:
        return "Input space points can not be used for parameter estimation";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::SeedingError e) {
  static SeedingErrorCategory c;
  return {static_cast<int>(e), c};
}
