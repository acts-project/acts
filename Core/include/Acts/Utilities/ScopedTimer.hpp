// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <chrono>
#include <iostream>

namespace Acts {

/// @brief A RAII timer class for measuring execution time of code blocks
///
/// ScopedTimer provides automatic timing of code blocks using RAII principles.
/// It starts timing when constructed and automatically logs the duration when
/// destroyed. This makes it ideal for measuring execution time of functions
/// or code blocks without manual start/stop calls.
///
/// Example usage:
/// @code
/// {
///   ScopedTimer timer("myFunction");
///   // ... code to measure ...
/// } // Timer automatically logs duration when block ends
/// @endcode
///
class ScopedTimer {
 public:
  /// @brief Construct a new Scoped Timer
  ///
  /// @param name Identifier for the timed block
  /// @param logger Logger instance to use for output
  /// @param lvl Logging level for the timing output
  explicit ScopedTimer(const std::string& name, const Logger& logger,
                       Logging::Level lvl = Logging::Level::INFO)
      : m_name(name), m_lvl(lvl), m_logger(&logger) {
    m_start = std::chrono::high_resolution_clock::now();
  }

  /// @brief Destructor that logs the execution time
  ///
  /// Automatically calculates and logs the duration between construction
  /// and destruction using the specified logger and level.
  ~ScopedTimer() {
    if (m_logger->doPrint(m_lvl)) {
      auto end = std::chrono::high_resolution_clock::now();
      auto duration =
          std::chrono::duration_cast<std::chrono::milliseconds>(end - m_start);
      std::ostringstream oss;
      oss << m_name << " took " << duration.count() << " ms";
      m_logger->log(m_lvl, oss.str());
    }
  }

 private:
  std::string m_name;      ///< Identifier for the timed block
  Logging::Level m_lvl;    ///< Logging level for output
  const Logger* m_logger;  ///< Logger instance for output
  std::chrono::high_resolution_clock::time_point
      m_start;  ///< Start time of the timer
};

}  // namespace Acts
