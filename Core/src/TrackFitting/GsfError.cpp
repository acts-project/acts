// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfError.hpp"

namespace {

class GsfErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "GsfError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::Experimental::GsfError;

    switch (static_cast<GsfError>(c)) {
      case GsfError::StartParametersNotOnStartSurface:
        return "Start parameters don't lie in the start surface";
      case GsfError::NoMeasurementStatesCreatedForward:
        return "No measurement states found in the forward pass";
      case GsfError::NoMeasurementStatesCreatedBackward:
        return "No measurement states found in the backward pass";
      case GsfError::NoMeasurementStatesCreatedFinal:
        return "No measurement states in the final trajectory";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::Experimental::make_error_code(
    Acts::Experimental::GsfError e) {
  static GsfErrorCategory c;
  return {static_cast<int>(e), c};
}
