// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Result.hpp"

#include <iostream>
#include <system_error>

// Example struct representing a track for demonstration purposes
struct Track {
  int id;
};

// Example struct representing track parameters
struct TrackParameters {
  double momentum;
};

//! [Error Enum Pattern]
enum class MyComponentError {
  ErrorValue1 = 1,  // All error values must be non-zero
  ErrorValue2,
  // ...
};

// Factory function to create std::error_code
std::error_code make_error_code(MyComponentError e);
//! [Error Enum Pattern]

//! [Usage with Result Type]
Acts::Result<Track> fitTrack(const TrackParameters& params) {
  if (params.momentum < 0) {
    return Acts::KalmanFitterError::UpdateFailed;
  }
  Track track{42};
  return track;  // Success
}

void processTrack() {
  TrackParameters params{100.0};

  // Usage
  auto result = fitTrack(params);
  if (!result.ok()) {
    std::error_code error = result.error();
    // Handle error
    std::cerr << "Error: " << error.message() << std::endl;
  } else {
    Track track = std::move(result).value();
    std::cout << "Track fitted successfully: " << track.id << std::endl;
  }
}
//! [Usage with Result Type]
