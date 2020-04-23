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
enum class CombinatorialKalmanFilterError {
  UpdateFailed = 1,
  SmoothFailed = 2,
  OutputConversionFailed = 3,
  SourcelinkSelectionFailed = 4,
  NoTracksFound = 5,
  PropagationReachesMaxSteps = 6
};

namespace detail {
// Define a custom error code category derived from std::error_category
class CombinatorialKalmanFilterErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final {
    return "CombinatorialKalmanFilterError";
  }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<CombinatorialKalmanFilterError>(c)) {
      case CombinatorialKalmanFilterError::UpdateFailed:
        return "Kalman update failed";
      case CombinatorialKalmanFilterError::SmoothFailed:
        return "Kalman smooth failed";
      case CombinatorialKalmanFilterError::OutputConversionFailed:
        return "Kalman output conversion failed";
      case CombinatorialKalmanFilterError::SourcelinkSelectionFailed:
        return "Source link selection failed";
      case CombinatorialKalmanFilterError::NoTracksFound:
        return "No track is found";
      case CombinatorialKalmanFilterError::PropagationReachesMaxSteps:
        return "Propagation reaches max steps before track finding is "
               "finished";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::CombinatorialKalmanFilterErrorCategory&
CombinatorialKalmanFilterErrorCategory() {
  static detail::CombinatorialKalmanFilterErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::CombinatorialKalmanFilterError e) {
  return {static_cast<int>(e), Acts::CombinatorialKalmanFilterErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::CombinatorialKalmanFilterError>
    : std::true_type {};
}  // namespace std
