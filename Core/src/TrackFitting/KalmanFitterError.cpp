// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <string>

namespace {

class KalmanFitterErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "KalmanFitterError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::KalmanFitterError;

    switch (static_cast<KalmanFitterError>(c)) {
      case KalmanFitterError::ForwardUpdateFailed:
        return "Kalman forward update failed";
      case KalmanFitterError::BackwardUpdateFailed:
        return "Kalman backward update failed";
      case KalmanFitterError::SmoothFailed:
        return "Kalman smooth failed";
      case KalmanFitterError::OutputConversionFailed:
        return "Kalman output conversion failed";
      case KalmanFitterError::NoMeasurementFound:
        return "No measurement detected during the propagation";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::KalmanFitterError e) {
  static KalmanFitterErrorCategory c;
  return {static_cast<int>(e), c};
}
