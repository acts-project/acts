// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_utils.hpp"
#include "detray/benchmarks/cpu/vector_benchmark.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

namespace detray {

namespace algebra {

template <detray::concepts::transform3D transform3_t>
void fill_random_trf(std::vector<transform3_t>&);

}  // namespace algebra

namespace benchmarks {

/// Benchmark for vector operations
template <detray::concepts::transform3D transform3_t>
struct transform3_bm : public vector_bm<typename transform3_t::vector3> {
 private:
  using base_type = vector_bm<typename transform3_t::vector3>;

 public:
  /// Prefix for the benchmark name
  static constexpr std::string_view bm_name{"TRANSFORM3"};

  std::vector<transform3_t> trfs;

  /// No default construction: Cannot prepare data
  transform3_bm() = delete;
  /// Construct from an externally provided configuration @param cfg
  explicit transform3_bm(benchmark_base::configuration cfg) : base_type{cfg} {
    trfs.reserve(static_cast<std::size_t>(this->m_cfg.n_samples()));

    detray::algebra::fill_random_trf(trfs);
  }
  transform3_bm(const transform3_bm& bm) = default;
  transform3_bm& operator=(transform3_bm& other) = default;

  /// Clear state
  ~transform3_bm() override { trfs.clear(); }

  std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{bm_name};
  }

  /// Benchmark case
  inline void operator()(::benchmark::State& state) const {
    using vector_t = typename transform3_t::vector3;
    using point_t = typename transform3_t::point3;

    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i = 0u; i < n_samples; ++i) {
        point_t result1 = this->trfs[i].point_to_global(this->a[i]);
        point_t result2 = this->trfs[i].point_to_local(this->a[i]);
        vector_t result3 = this->trfs[i].vector_to_global(this->a[i]);
        vector_t result4 = this->trfs[i].vector_to_local(this->a[i]);

        ::benchmark::DoNotOptimize(result1);
        ::benchmark::DoNotOptimize(result2);
        ::benchmark::DoNotOptimize(result3);
        ::benchmark::DoNotOptimize(result4);
      }
    }
  }
};

// Macro for defining all transform benchmark types
#define DETRAY_DEFINE_TRANSFORM_BENCH(PLUGIN)               \
  using trf_f_t = transform3_bm<PLUGIN::transform3<float>>; \
  using trf_d_t = transform3_bm<PLUGIN::transform3<double>>;

// Macro for registering all transform benchmarks
#define DETRAY_REGISTER_TRANSFORM_BENCH(CFGS, CFGD)                 \
  detray::benchmarks::register_benchmark<trf_f_t>(CFGS, "_SINGLE"); \
  detray::benchmarks::register_benchmark<trf_d_t>(CFGD, "_DOUBLE");

}  // namespace benchmarks

}  // namespace detray
