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
  using clock_type = std::chrono::high_resolution_clock;

  /// @brief Construct a new Scoped Timer
  ///
  /// @param name Identifier for the timed block
  /// @param logger Logger instance to use for output
  /// @param lvl Logging level for the timing output
  explicit ScopedTimer(const std::string& name, const Logger& logger,
                       Logging::Level lvl = Logging::Level::INFO);

  /// @brief Destructor that logs the execution time
  ///
  /// Automatically calculates and logs the duration between construction
  /// and destruction using the specified logger and level.
  ~ScopedTimer();

  ScopedTimer(const ScopedTimer&) = delete;
  ScopedTimer& operator=(const ScopedTimer&) = delete;
  ScopedTimer(ScopedTimer&&) = delete;
  ScopedTimer& operator=(ScopedTimer&&) = delete;

 private:
  std::string m_name;              ///< Identifier for the timed block
  Logging::Level m_lvl;            ///< Logging level for output
  const Logger* m_logger;          ///< Logger instance for output
  clock_type::time_point m_start;  ///< Start time of the timer
};

/// @brief A timer class that measures and averages execution times of multiple samples
///
/// This class provides functionality to measure execution times of code blocks
/// and calculate statistics (mean, standard deviation) across multiple samples.
/// It uses RAII through the Sample class to automatically record timing
/// information.
class AveragingScopedTimer {
 public:
  using clock_type = std::chrono::high_resolution_clock;

  /// @brief RAII wrapper class for measuring individual timing samples
  ///
  /// When constructed, starts a timer. When destroyed, automatically records
  /// the duration to the parent AveragingScopedTimer.
  class Sample {
   public:
    /// @brief Construct a new sample and start timing
    explicit Sample(AveragingScopedTimer& parent);
    /// @brief Record the duration when destroyed
    ~Sample();
    Sample(const Sample&) = delete;
    Sample& operator=(const Sample&) = delete;
    /// @brief Move constructor that transfers ownership of timing to new sample
    Sample(Sample&& /*other*/);
    Sample& operator=(Sample&&) = delete;

   private:
    AveragingScopedTimer* m_parent;
    clock_type::time_point m_start;
  };

  /// @brief Construct a new AveragingScopedTimer
  ///
  /// @param name Name of the timer for logging
  /// @param logger Logger instance to use for output
  /// @param lvl Logging level for timing output
  explicit AveragingScopedTimer(const std::string& name, const Logger& logger,
                                Logging::Level lvl = Logging::Level::INFO);

  /// @brief Destroy the AveragingScopedTimer and log statistics
  ///
  /// Outputs total duration and per-sample statistics (mean ± stddev) if
  /// logging is enabled at the configured level.
  ~AveragingScopedTimer();
  AveragingScopedTimer(const AveragingScopedTimer&) = delete;
  AveragingScopedTimer& operator=(const AveragingScopedTimer&) = delete;
  AveragingScopedTimer(AveragingScopedTimer&&) = delete;
  AveragingScopedTimer& operator=(AveragingScopedTimer&&) = delete;

  /// @brief Create a new timing sample
  ///
  /// @return Sample RAII wrapper for measuring a single timing sample
  Sample sample();

  friend class Sample;

 private:
  /// @brief Add a timing sample to the statistics
  ///
  /// @param duration Duration of the sample in nanoseconds
  void addSample(std::chrono::nanoseconds duration);

  double m_sumDuration = 0;  ///< Sum of all sample durations
  double m_sumDurationSquared =
      0;  ///< Sum of squared durations for stddev calculation
  std::size_t m_nSamples = 0;  ///< Number of samples recorded

  std::string m_name;      ///< Name of the timer for logging
  Logging::Level m_lvl;    ///< Logging level for output
  const Logger* m_logger;  ///< Logger instance for output
};
}  // namespace Acts
