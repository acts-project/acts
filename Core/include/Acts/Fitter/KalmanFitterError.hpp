// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace Acts {
// This is the custom error code enum
enum class KalmanFitterError {
  ForwardUpdateFailed = 1,
  BackwardUpdateFailed = 2,
  SmoothFailed = 3,
  OutputConversionFailed = 4,
  NoMeasurementFound = 5
};

namespace detail {
// Define a custom error code category derived from std::error_category
class KalmanFitterErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "KalmanFitterError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
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
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::KalmanFitterErrorCategory&
KalmanFitterErrorCategory() {
  static detail::KalmanFitterErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::KalmanFitterError e) {
  return {static_cast<int>(e), Acts::KalmanFitterErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::KalmanFitterError> : std::true_type {};
}  // namespace std
