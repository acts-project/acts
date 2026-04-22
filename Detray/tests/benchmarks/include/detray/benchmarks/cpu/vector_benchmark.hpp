// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/benchmark_utils.hpp"

// System include(s)
#include <chrono>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

namespace detray {

namespace algebra {

template <concepts::vector vector_t>
void fill_random_vec(std::vector<vector_t> &);

}  // namespace algebra

namespace benchmarks {

/// Benchmark for vector operations
template <concepts::vector vector_t>
struct vector_bm : public benchmark_base {
  /// Prefix for the benchmark name
  static constexpr std::string_view name{"BM_VECTOR"};

  std::vector<vector_t> a;
  std::vector<vector_t> b;
  std::vector<vector_t> results;

  /// No default construction: Cannot prepare data
  vector_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  explicit vector_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {
    const auto n_data{static_cast<std::size_t>(this->m_cfg.n_samples())};

    a.reserve(n_data);
    b.reserve(n_data);

    detray::algebra::fill_random_vec(a);
    detray::algebra::fill_random_vec(b);
  }
  vector_bm(const vector_bm &bm) = default;
  vector_bm &operator=(vector_bm &other) = default;

  /// Clear state
  ~vector_bm() override {
    a.clear();
    b.clear();
  }
};

/// Benchmark unary operations on vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename unaryOP>
  requires std::invocable<unaryOP, vector_t<scalar_t>>
struct vector_unaryOP_bm : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_unaryOP_bm() = delete;
  explicit vector_unaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  vector_unaryOP_bm(const vector_unaryOP_bm &bm) = default;
  vector_unaryOP_bm &operator=(vector_unaryOP_bm &other) = default;

  std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{unaryOP::name};
  }

  inline void operator()(::benchmark::State &state) const {
    using result_t = std::invoke_result_t<unaryOP, vector_t<scalar_t>>;

    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = unaryOP{}(this->a[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

/// Benchmark binary operations on vectors
template <template <typename> class vector_t, concepts::scalar scalar_t,
          typename binaryOP>
  requires std::invocable<binaryOP, vector_t<scalar_t>, vector_t<scalar_t>>
struct vector_binaryOP_bm : public vector_bm<vector_t<scalar_t>> {
  using base_type = vector_bm<vector_t<scalar_t>>;

  vector_binaryOP_bm() = delete;
  explicit vector_binaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  vector_binaryOP_bm(const vector_binaryOP_bm &bm) = default;
  vector_binaryOP_bm &operator=(vector_binaryOP_bm &other) = default;

  std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{binaryOP::name};
  }

  inline void operator()(::benchmark::State &state) const {
    using result_t =
        std::invoke_result_t<binaryOP, vector_t<scalar_t>, vector_t<scalar_t>>;

    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i{0}; i < n_samples; ++i) {
        result_t result = binaryOP{}(this->a[i], this->b[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

// Operations to be benchmarked
namespace bench_op::vector {

struct add {
  static constexpr std::string_view name{"ADD"};
  template <concepts::vector vector_t>
  constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
    return a + b;
  }
};
struct sub {
  static constexpr std::string_view name{"SUB"};
  template <concepts::vector vector_t>
  constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
    return a - b;
  }
};
struct dot {
  static constexpr std::string_view name{"DOT"};
  template <concepts::vector vector_t>
  constexpr detray::traits::scalar_t<vector_t> operator()(
      const vector_t &a, const vector_t &b) const {
    return detray::vector::dot(a, b);
  }
};
struct cross {
  static constexpr std::string_view name{"CROSS"};
  template <concepts::vector vector_t>
  constexpr vector_t operator()(const vector_t &a, const vector_t &b) const {
    return detray::vector::cross(a, b);
  }
};

// Macro for declaring vector unary ops
#define DETRAY_BENCH_VECTOR(OP, NAME, RES)              \
  struct OP {                                           \
    static constexpr std::string_view name{#NAME};      \
    template <concepts::vector vector_t>                \
    constexpr RES operator()(const vector_t &a) const { \
      return detray::vector::OP(a);                     \
    }                                                   \
  };

DETRAY_BENCH_VECTOR(phi, PHI, detray::traits::scalar_t<vector_t>)
DETRAY_BENCH_VECTOR(theta, THETA, detray::traits::scalar_t<vector_t>)
DETRAY_BENCH_VECTOR(eta, ETA, detray::traits::scalar_t<vector_t>)
DETRAY_BENCH_VECTOR(perp, PERP, detray::traits::scalar_t<vector_t>)
DETRAY_BENCH_VECTOR(norm, NORM, detray::traits::scalar_t<vector_t>)
DETRAY_BENCH_VECTOR(normalize, NORMALIZE, vector_t)

}  // namespace bench_op::vector

// Macro for defining all vector benchmark types
#define DETRAY_DEFINE_VECTOR_BENCH(PLUGIN)                                    \
  using phi_f_t =                                                             \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::phi>;       \
  using theta_f_t =                                                           \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::theta>;     \
  using perp_f_t =                                                            \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::perp>;      \
  using norm_f_t =                                                            \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::norm>;      \
  using eta_f_t =                                                             \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::eta>;       \
                                                                              \
  using add_f_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::add>;      \
  using sub_f_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::sub>;      \
  using dot_f_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::dot>;      \
  using cross_f_t =                                                           \
      vector_binaryOP_bm<PLUGIN::vector3, float, bench_op::vector::cross>;    \
  using normlz_f_t =                                                          \
      vector_unaryOP_bm<PLUGIN::vector3, float, bench_op::vector::normalize>; \
                                                                              \
  using phi_d_t =                                                             \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::phi>;      \
  using theta_d_t =                                                           \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::theta>;    \
  using perp_d_t =                                                            \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::perp>;     \
  using norm_d_t =                                                            \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::norm>;     \
  using eta_d_t =                                                             \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::eta>;      \
                                                                              \
  using add_d_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::add>;     \
  using sub_d_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::sub>;     \
  using dot_d_t =                                                             \
      vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::dot>;     \
  using cross_d_t =                                                           \
      vector_binaryOP_bm<PLUGIN::vector3, double, bench_op::vector::cross>;   \
  using normlz_d_t =                                                          \
      vector_unaryOP_bm<PLUGIN::vector3, double, bench_op::vector::normalize>;

// Macro for registering all vector benchmarks
#define DETRAY_REGISTER_VECTOR_BENCH(CFGS, CFGD)                       \
  detray::benchmarks::register_benchmark<add_f_t>(CFGS, "_SINGLE");    \
  detray::benchmarks::register_benchmark<add_d_t>(CFGD, "_DOUBLE");    \
  detray::benchmarks::register_benchmark<sub_f_t>(CFGS, "_SINGLE");    \
  detray::benchmarks::register_benchmark<sub_d_t>(CFGD, "_DOUBLE");    \
  detray::benchmarks::register_benchmark<dot_f_t>(CFGS, "_SINGLE");    \
  detray::benchmarks::register_benchmark<dot_d_t>(CFGD, "_DOUBLE");    \
  detray::benchmarks::register_benchmark<cross_f_t>(CFGS, "_SINGLE");  \
  detray::benchmarks::register_benchmark<cross_d_t>(CFGD, "_DOUBLE");  \
  detray::benchmarks::register_benchmark<normlz_f_t>(CFGS, "_SINGLE"); \
  detray::benchmarks::register_benchmark<normlz_d_t>(CFGD, "_DOUBLE"); \
                                                                       \
  detray::benchmarks::register_benchmark<phi_f_t>(CFGS, "_SINGLE");    \
  detray::benchmarks::register_benchmark<phi_d_t>(CFGD, "_DOUBLE");    \
  detray::benchmarks::register_benchmark<theta_f_t>(CFGS, "_SINGLE");  \
  detray::benchmarks::register_benchmark<theta_d_t>(CFGD, "_DOUBLE");  \
  detray::benchmarks::register_benchmark<perp_f_t>(CFGS, "_SINGLE");   \
  detray::benchmarks::register_benchmark<perp_d_t>(CFGD, "_DOUBLE");   \
  detray::benchmarks::register_benchmark<norm_f_t>(CFGS, "_SINGLE");   \
  detray::benchmarks::register_benchmark<norm_d_t>(CFGD, "_DOUBLE");   \
  detray::benchmarks::register_benchmark<eta_f_t>(CFGS, "_SINGLE");    \
  detray::benchmarks::register_benchmark<eta_d_t>(CFGD, "_DOUBLE");

}  // namespace benchmarks

}  // namespace detray
