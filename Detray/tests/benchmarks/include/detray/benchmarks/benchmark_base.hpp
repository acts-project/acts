// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <ostream>
#include <string>

namespace detray::benchmarks {

/// Benchmark configuration type
struct configuration {
  /// Prefix for the benchmark name
  std::string m_name{"DETRAY_BM_"};
  /// Size of data sample to be used in benchmark
  int m_samples{100};
  /// Run a number of operations before the benchmark
  bool m_warmup = true;
  // Size of data in warm-up round
  int m_n_warmup{static_cast<int>(0.1 * static_cast<double>(m_samples))};
  // Sleep after building data sample
  bool m_sleep = false;
  // Number of seconds to sleep
  std::size_t m_n_sleep{1u};

  /// Setters
  /// @{
  configuration& name(const std::string_view n) {
    m_name = n;
    return *this;
  }
  configuration& n_samples(int n) {
    m_samples = n;
    return *this;
  }
  configuration& do_warmup(bool b) {
    m_warmup = b;
    return *this;
  }
  configuration& n_warmup(int n) {
    m_n_warmup = n;
    m_warmup = true;
    return *this;
  }
  configuration& do_sleep(bool b) {
    m_sleep = b;
    return *this;
  }
  configuration& n_sleep(std::size_t n) {
    m_n_sleep = n;
    m_sleep = true;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  const std::string& name() const { return m_name; }
  constexpr int n_samples() const { return m_samples; }
  constexpr bool do_warmup() const { return m_warmup; }
  constexpr int n_warmup() const { return m_n_warmup; }
  constexpr bool do_sleep() const { return m_sleep; }
  constexpr std::size_t n_sleep() const { return m_n_sleep; }
  /// @}

 private:
  /// Print the benchmark setup
  friend std::ostream& operator<<(std::ostream& os, const configuration& cfg) {
    os << " -> running:\t " << cfg.n_samples() << " samples" << std::endl;
    if (cfg.do_warmup()) {
      os << " -> warmup: \t " << cfg.n_warmup() << " samples" << std::endl;
    }
    if (cfg.do_sleep()) {
      os << " -> cool down:\t " << cfg.n_sleep() << "s" << std::endl;
    }
    os << std::endl;
    return os;
  }
};

/// Base type for detray benchmarks using google benchmark
struct benchmark_base {
  using configuration = detray::benchmarks::configuration;

  /// Default construction
  benchmark_base() = default;

  /// Construct from an externally provided configuration @param cfg
  explicit benchmark_base(const configuration& cfg) : m_cfg{cfg} {}

  /// @returns the benchmark configuration
  const configuration& config() const { return m_cfg; }
  configuration& config() { return m_cfg; }

  /// Default destructor
  virtual ~benchmark_base() = default;

  /// @returns the benchmark name
  virtual std::string name() const { return m_cfg.name(); };

 protected:
  /// The benchmark configuration
  configuration m_cfg{};
};

}  // namespace detray::benchmarks
