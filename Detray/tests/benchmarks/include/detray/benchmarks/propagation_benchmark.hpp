// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/propagator/propagation_config.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"

// System include(s)
#include <string>
#include <string_view>

namespace detray::benchmarks {

/// Configuration for propagation benchmarks
struct propagation_benchmark_config : public detray::benchmarks::configuration {
  /// Propagation configuration
  propagation::config m_propagation{};

  /// Default construction
  propagation_benchmark_config() = default;

  /// Construct from a base configuration
  explicit propagation_benchmark_config(
      const detray::benchmarks::configuration& bench_cfg)
      : detray::benchmarks::configuration(bench_cfg) {}

  /// Construct from a base configuration
  explicit propagation_benchmark_config(
      const detray::benchmarks::configuration& bench_cfg,
      const propagation::config& prop_cfg)
      : detray::benchmarks::configuration(bench_cfg), m_propagation{prop_cfg} {}

  /// Getters
  /// @{
  const propagation::config& propagation() const { return m_propagation; }
  propagation::config& propagation() { return m_propagation; }
  /// @}
};

/// Base type for propagation benchmarks
template <concepts::algebra algebra_t>
struct propagation_benchmark : public benchmark_base {
  /// Detector dependent types
  using scalar_t = dscalar<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

  /// Local configuration type
  using configuration = propagation_benchmark_config;

  /// Default construction
  propagation_benchmark() = default;

  /// Construct from an externally provided configuration @param cfg
  explicit propagation_benchmark(const configuration& cfg)
      : benchmark_base(cfg), m_prop_cfg{cfg.propagation()} {}

  /// @return the benchmark configuration
  const propagation::config& propagation() const { return m_prop_cfg; }
  propagation::config& propagation() { return m_prop_cfg; }

 protected:
  /// The benchmark configuration
  propagation::config m_prop_cfg{};
};

}  // namespace detray::benchmarks
