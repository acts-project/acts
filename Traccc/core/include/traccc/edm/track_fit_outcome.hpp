// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <cstdint>

namespace traccc {

/// Possible outcomes of a track fit
enum class track_fit_outcome : std::uint16_t {
  UNKNOWN,
  SUCCESS,
  FAILURE_NON_POSITIVE_NDF,
  FAILURE_NOT_ALL_FITTED,
  FAILURE_NOT_ALL_SMOOTHED,
  FAILURE_FITTER,
  FAILURE_SMOOTHER,
  FAILURE_FORWARD_PROPAGATION,
  FAILURE_BACKWARD_PROPAGATION,
  MAX_OUTCOME
};

/// Convert a fit outcome code into a human readable string
struct fit_outcome_debug_msg {
  TRACCC_HOST std::string operator()() const {
    switch (m_code) {
      using enum track_fit_outcome;
      case UNKNOWN: {
        return "Unknown";
      }
      case SUCCESS: {
        return "Success";
      }
      case FAILURE_NON_POSITIVE_NDF: {
        return "Negative NDF";
      }
      case FAILURE_NOT_ALL_FITTED: {
        return "Not all track states fitted";
      }
      case FAILURE_NOT_ALL_SMOOTHED: {
        return "Not all track states smoothed";
      }
      case FAILURE_FITTER: {
        return "Fitter failure";
      }
      case FAILURE_SMOOTHER: {
        return "Smoother failure";
      }
      case FAILURE_FORWARD_PROPAGATION: {
        return "Forward propagation failure";
      }
      case FAILURE_BACKWARD_PROPAGATION: {
        return "Backward propagation failure";
      }
      default: {
        return "Not a defined fit outcome code: " +
               std::to_string(static_cast<int>(m_code));
      }
    }
  }

  track_fit_outcome m_code{track_fit_outcome::MAX_OUTCOME};
};

}  // namespace traccc
