// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/benchmark_utils.hpp"

// System include(s)
#include <string_view>
#include <vector>

namespace detray {

namespace algebra {

template <detray::concepts::matrix matrix_t>
void fill_random_matrix(std::vector<matrix_t>&);

template <detray::concepts::vector vector_t>
void fill_random_vec(std::vector<vector_t>&);

}  // namespace algebra

namespace benchmarks {

/// Benchmark for matrix operations
template <detray::concepts::matrix matrix_t>
struct matrix_bm : public benchmark_base {
  /// Prefix for the benchmark name
  static constexpr std::string_view name{"BM_MATRIX"};

  std::vector<matrix_t> a;
  std::vector<matrix_t> b;

  /// No default construction: Cannot prepare data
  matrix_bm() = delete;

  /// Construct from an externally provided configuration @param cfg
  explicit matrix_bm(benchmark_base::configuration cfg) : benchmark_base{cfg} {
    const auto n_data{static_cast<std::size_t>(this->m_cfg.n_samples())};

    a.reserve(n_data);
    b.reserve(n_data);

    detray::algebra::fill_random_matrix(a);
    detray::algebra::fill_random_matrix(b);
  }

  matrix_bm(const matrix_bm& bm) = default;
  matrix_bm& operator=(matrix_bm& other) = default;

  /// Clear state
  ~matrix_bm() override {
    a.clear();
    b.clear();
  }
};

/// Benchmark operations on a single matrix (transpose, inverse etc)
template <detray::concepts::matrix matrix_t, typename unaryOP>
  requires std::invocable<unaryOP, matrix_t>
struct matrix_unaryOP_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_unaryOP_bm() = delete;
  explicit matrix_unaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  matrix_unaryOP_bm(const matrix_unaryOP_bm& bm) = default;
  matrix_unaryOP_bm& operator=(matrix_unaryOP_bm& other) = default;

  std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{unaryOP::name};
  }

  inline void operator()(::benchmark::State& state) const {
    using result_t = std::invoke_result_t<unaryOP, matrix_t>;

    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i = 0u; i < n_samples; ++i) {
        result_t result = unaryOP{}(this->a[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

/// Benchmark elementwise addition/subtraction/multiplication of matrices
template <detray::concepts::matrix matrix_t, typename binaryOP>
  requires std::invocable<binaryOP, matrix_t, matrix_t>
struct matrix_binaryOP_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  matrix_binaryOP_bm() = delete;
  explicit matrix_binaryOP_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {}
  matrix_binaryOP_bm(const matrix_binaryOP_bm& bm) = default;
  matrix_binaryOP_bm& operator=(matrix_binaryOP_bm& other) = default;

  std::string name() const override {
    return std::string{base_type::name} + "_" + std::string{binaryOP::name};
  }

  inline void operator()(::benchmark::State& state) const {
    using result_t = std::invoke_result_t<binaryOP, matrix_t, matrix_t>;

    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i = 0u; i < n_samples; ++i) {
        result_t result = binaryOP{}(this->a[i], this->b[i]);
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

/// Benchmark matrix vector multiplication
template <detray::concepts::matrix matrix_t, concepts::vector vector_t>
struct matrix_vector_bm : public matrix_bm<matrix_t> {
  using base_type = matrix_bm<matrix_t>;

  std::vector<vector_t> v;

  matrix_vector_bm() = delete;
  explicit matrix_vector_bm(benchmark_base::configuration cfg)
      : base_type{cfg} {
    v.reserve(static_cast<std::size_t>(this->m_cfg.n_samples()));

    detray::algebra::fill_random_vec(v);
  }
  matrix_vector_bm(const matrix_vector_bm& bm) = default;
  matrix_vector_bm& operator=(matrix_vector_bm& other) = default;

  /// Clear state
  ~matrix_vector_bm() override { v.clear(); }

  std::string name() const override {
    return std::string{base_type::name} + "_VECTOR";
  }

  inline void operator()(::benchmark::State& state) const {
    const auto n_samples{static_cast<std::size_t>(this->m_cfg.n_samples())};

    // Run the benchmark
    for (auto _ : state) {
      for (std::size_t i = 0u; i < n_samples; ++i) {
        vector_t result = this->a[i] * this->v[i];
        ::benchmark::DoNotOptimize(result);
      }
    }
  }
};

// Functions to be benchmarked
namespace bench_op::matrix {

struct add {
  static constexpr std::string_view name{"ADD"};
  template <detray::concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a + b;
  }
};
struct sub {
  static constexpr std::string_view name{"SUB"};
  template <detray::concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a - b;
  }
};
struct mul {
  static constexpr std::string_view name{"MUL"};
  template <detray::concepts::matrix matrix_t>
  constexpr matrix_t operator()(const matrix_t& a, const matrix_t& b) const {
    return a * b;
  }
};
struct transpose {
  static constexpr std::string_view name{"TRANSPOSE"};
  template <detray::concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return detray::matrix::transpose(a);
  }
};
struct determinant {
  static constexpr std::string_view name{"DETERMINANT"};
  template <detray::concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return detray::matrix::determinant(a);
  }
};
struct invert {
  static constexpr std::string_view name{"INVERT"};
  template <detray::concepts::matrix matrix_t>
  constexpr auto operator()(const matrix_t& a) const {
    return detray::matrix::inverse(a);
  }
};

}  // namespace bench_op::matrix

// Macro for defining all matrix benchmark types
#define DETRAY_DEFINE_MATRIX_BENCH(PLUGIN)                                     \
  using mat44_transp_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 4, 4>, \
                                             bench_op::matrix::transpose>;     \
  using mat44_transp_d_t =                                                     \
      matrix_unaryOP_bm<PLUGIN::matrix_type<double, 4, 4>,                     \
                        bench_op::matrix::transpose>;                          \
  using mat66_transp_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 6, 6>, \
                                             bench_op::matrix::transpose>;     \
  using mat66_transp_d_t =                                                     \
      matrix_unaryOP_bm<PLUGIN::matrix_type<double, 6, 6>,                     \
                        bench_op::matrix::transpose>;                          \
  using mat88_transp_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 8, 8>, \
                                             bench_op::matrix::transpose>;     \
  using mat88_transp_d_t =                                                     \
      matrix_unaryOP_bm<PLUGIN::matrix_type<double, 8, 8>,                     \
                        bench_op::matrix::transpose>;                          \
                                                                               \
  using mat44_inv_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 4, 4>,    \
                                          bench_op::matrix::invert>;           \
  using mat44_inv_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 4, 4>,   \
                                          bench_op::matrix::invert>;           \
  using mat66_inv_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 6, 6>,    \
                                          bench_op::matrix::invert>;           \
  using mat66_inv_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 6, 6>,   \
                                          bench_op::matrix::invert>;           \
  using mat88_inv_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 8, 8>,    \
                                          bench_op::matrix::invert>;           \
  using mat88_inv_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 8, 8>,   \
                                          bench_op::matrix::invert>;           \
                                                                               \
  using mat44_det_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 4, 4>,    \
                                          bench_op::matrix::determinant>;      \
  using mat44_det_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 4, 4>,   \
                                          bench_op::matrix::determinant>;      \
  using mat66_det_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 6, 6>,    \
                                          bench_op::matrix::determinant>;      \
  using mat66_det_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 6, 6>,   \
                                          bench_op::matrix::determinant>;      \
  using mat88_det_f_t = matrix_unaryOP_bm<PLUGIN::matrix_type<float, 8, 8>,    \
                                          bench_op::matrix::determinant>;      \
  using mat88_det_d_t = matrix_unaryOP_bm<PLUGIN::matrix_type<double, 8, 8>,   \
                                          bench_op::matrix::determinant>;      \
                                                                               \
  using mat44_add_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 4, 4>,   \
                                           bench_op::matrix::add>;             \
  using mat44_add_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 4, 4>,  \
                                           bench_op::matrix::add>;             \
  using mat66_add_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 6, 6>,   \
                                           bench_op::matrix::add>;             \
  using mat66_add_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 6, 6>,  \
                                           bench_op::matrix::add>;             \
  using mat88_add_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 8, 8>,   \
                                           bench_op::matrix::add>;             \
  using mat88_add_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 8, 8>,  \
                                           bench_op::matrix::add>;             \
                                                                               \
  using mat44_mul_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 4, 4>,   \
                                           bench_op::matrix::mul>;             \
  using mat44_mul_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 4, 4>,  \
                                           bench_op::matrix::mul>;             \
  using mat66_mul_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 6, 6>,   \
                                           bench_op::matrix::mul>;             \
  using mat66_mul_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 6, 6>,  \
                                           bench_op::matrix::mul>;             \
  using mat88_mul_f_t = matrix_binaryOP_bm<PLUGIN::matrix_type<float, 8, 8>,   \
                                           bench_op::matrix::mul>;             \
  using mat88_mul_d_t = matrix_binaryOP_bm<PLUGIN::matrix_type<double, 8, 8>,  \
                                           bench_op::matrix::mul>;             \
                                                                               \
  using mat44_vec_f_t = matrix_vector_bm<PLUGIN::matrix_type<float, 4, 4>,     \
                                         PLUGIN::vector_type<float, 4>>;       \
  using mat44_vec_d_t = matrix_vector_bm<PLUGIN::matrix_type<double, 4, 4>,    \
                                         PLUGIN::vector_type<double, 4>>;      \
  using mat66_vec_f_t = matrix_vector_bm<PLUGIN::matrix_type<float, 6, 6>,     \
                                         PLUGIN::vector_type<float, 6>>;       \
  using mat66_vec_d_t = matrix_vector_bm<PLUGIN::matrix_type<double, 6, 6>,    \
                                         PLUGIN::vector_type<double, 6>>;      \
  using mat88_vec_f_t = matrix_vector_bm<PLUGIN::matrix_type<float, 8, 8>,     \
                                         PLUGIN::vector_type<float, 8>>;       \
  using mat88_vec_d_t = matrix_vector_bm<PLUGIN::matrix_type<double, 8, 8>,    \
                                         PLUGIN::vector_type<double, 8>>;

// Macro for registering all matrix benchmarks
#define DETRAY_REGISTER_MATRIX_BENCH(CFGS, CFGD)                              \
  detray::benchmarks::register_benchmark<mat44_transp_f_t>(CFGS,              \
                                                           "_4x4_SINGLE");    \
  detray::benchmarks::register_benchmark<mat44_transp_d_t>(CFGD,              \
                                                           "_4x4_DOUBLE");    \
  detray::benchmarks::register_benchmark<mat66_transp_f_t>(CFGS,              \
                                                           "_6x6_SINGLE");    \
  detray::benchmarks::register_benchmark<mat66_transp_d_t>(CFGD,              \
                                                           "_6x6_DOUBLE");    \
  detray::benchmarks::register_benchmark<mat88_transp_f_t>(CFGS,              \
                                                           "_8x8_SINGLE");    \
  detray::benchmarks::register_benchmark<mat88_transp_d_t>(CFGD,              \
                                                           "_8x8_DOUBLE");    \
                                                                              \
  detray::benchmarks::register_benchmark<mat44_inv_f_t>(CFGS, "_4x4_SINGLE"); \
  detray::benchmarks::register_benchmark<mat44_inv_d_t>(CFGD, "_4x4_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat66_inv_f_t>(CFGS, "_6x6_SINGLE"); \
  detray::benchmarks::register_benchmark<mat66_inv_d_t>(CFGD, "_6x6_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat88_inv_f_t>(CFGS, "_8x8_SINGLE"); \
  detray::benchmarks::register_benchmark<mat88_inv_d_t>(CFGD, "_8x8_DOUBLE"); \
                                                                              \
  detray::benchmarks::register_benchmark<mat44_det_f_t>(CFGS, "_4x4_SINGLE"); \
  detray::benchmarks::register_benchmark<mat44_det_d_t>(CFGD, "_4x4_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat66_det_f_t>(CFGS, "_6x6_SINGLE"); \
  detray::benchmarks::register_benchmark<mat66_det_d_t>(CFGD, "_6x6_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat88_det_f_t>(CFGS, "_8x8_SINGLE"); \
  detray::benchmarks::register_benchmark<mat88_det_d_t>(CFGD, "_8x8_DOUBLE"); \
                                                                              \
  detray::benchmarks::register_benchmark<mat44_add_f_t>(CFGS, "_4x4_SINGLE"); \
  detray::benchmarks::register_benchmark<mat44_add_d_t>(CFGD, "_4x4_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat66_add_f_t>(CFGS, "_6x6_SINGLE"); \
  detray::benchmarks::register_benchmark<mat66_add_d_t>(CFGD, "_6x6_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat88_add_f_t>(CFGS, "_8x8_SINGLE"); \
  detray::benchmarks::register_benchmark<mat88_add_d_t>(CFGD, "_8x8_DOUBLE"); \
                                                                              \
  detray::benchmarks::register_benchmark<mat44_mul_f_t>(CFGS, "_4x4_SINGLE"); \
  detray::benchmarks::register_benchmark<mat44_mul_d_t>(CFGD, "_4x4_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat66_mul_f_t>(CFGS, "_6x6_SINGLE"); \
  detray::benchmarks::register_benchmark<mat66_mul_d_t>(CFGD, "_6x6_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat88_mul_f_t>(CFGS, "_8x8_SINGLE"); \
  detray::benchmarks::register_benchmark<mat88_mul_d_t>(CFGD, "_8x8_DOUBLE"); \
                                                                              \
  detray::benchmarks::register_benchmark<mat44_vec_f_t>(CFGS, "_4x4_SINGLE"); \
  detray::benchmarks::register_benchmark<mat44_vec_d_t>(CFGD, "_4x4_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat66_vec_f_t>(CFGS, "_6x6_SINGLE"); \
  detray::benchmarks::register_benchmark<mat66_vec_d_t>(CFGD, "_6x6_DOUBLE"); \
  detray::benchmarks::register_benchmark<mat88_vec_f_t>(CFGS, "_8x8_SINGLE"); \
  detray::benchmarks::register_benchmark<mat88_vec_d_t>(CFGD,                 \
                                                        "_8x8_"               \
                                                        "double");

}  // namespace benchmarks

}  // namespace detray
