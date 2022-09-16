// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <chrono>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

namespace Acts {

/// @struct ExaTrkXTime
///
/// @brief Collection of timing information of the Exa.TrkX algorithm
struct ExaTrkXTime {
  float embedding = 0.0;
  float building = 0.0;
  float filtering = 0.0;
  float gnn = 0.0;
  float labeling = 0.0;
  float total = 0.0;

  void summary(LoggerWrapper& logger) const {
    ACTS_VERBOSE("1) embedding: " << embedding);
    ACTS_VERBOSE("2) building: " << building);
    ACTS_VERBOSE("3) filtering: " << filtering);
    ACTS_VERBOSE("4) gnn: " << gnn);
    ACTS_VERBOSE("5) labeling: " << labeling);
    ACTS_VERBOSE("6) total: " << total);
  }
};

/// @class ExaTrkXTimer
///
/// A timer to allow easy timing
class ExaTrkXTimer {
 public:
  ExaTrkXTimer(bool disabled = false) : m_disabled(disabled) {}

  void start() {
    if (not m_disabled) {
      m_start = std::chrono::high_resolution_clock::now();
      m_running = true;
    }
  }
  void stop() {
    if (not m_disabled) {
      m_end = std::chrono::high_resolution_clock::now();
      m_running = false;
    }
  }
  double stopAndGetElapsedTime() {
    stop();
    return elapsedSeconds();
  }
  double elapsed() {
    if (not m_disabled) {
      std::chrono::time_point<std::chrono::high_resolution_clock> end;
      if (m_running) {
        end = std::chrono::high_resolution_clock::now();
      } else {
        end = m_end;
      }
      return std::chrono::duration<double, std::milli>(end - m_start).count();
    } else {
      return 0.0;
    }
  }
  double elapsedSeconds() { return elapsed() / 1000.0; }

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
  bool m_running = false;
  bool m_disabled;
};

}  // namespace Acts
