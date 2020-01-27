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
enum class TrackFinderError {
  UpdateFailed = 1,
  SmoothFailed = 2,
  OutputConversionFailed = 3,
  SourcelinkSelectionFailed = 4,
  NoTracksFound = 5,
  PropagationReachesLimit = 6
};

namespace detail {
// Define a custom error code category derived from std::error_category
class TrackFinderErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "TrackFinderError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<TrackFinderError>(c)) {
      case TrackFinderError::UpdateFailed:
        return "Kalman update failed";
      case TrackFinderError::SmoothFailed:
        return "Kalman smooth failed";
      case TrackFinderError::OutputConversionFailed:
        return "Kalman output conversion failed";
      case TrackFinderError::SourcelinkSelectionFailed:
        return "Source link selection failed";
      case TrackFinderError::NoTracksFound:
        return "No track is found";
      case TrackFinderError::PropagationReachesLimit:
        return "Propagation reaches path limit before track finding is "
               "finished";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::TrackFinderErrorCategory&
TrackFinderErrorCategory() {
  static detail::TrackFinderErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::TrackFinderError e) {
  return {static_cast<int>(e), Acts::TrackFinderErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::TrackFinderError> : std::true_type {};
}  // namespace std
