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

// Project include(s)
// clang-format off
#include "algebra/vc_soa.hpp"
#include "algebra/utils/data_generator.hpp"
#include "detray/benchmarks/cpu/matrix_benchmark.hpp"
#include "detray/benchmarks/cpu/transform_benchmark.hpp"
#include "detray/benchmarks/cpu/vector_benchmark.hpp"
// clang-format on
#include "detray/benchmarks/benchmark_utils.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>

using namespace detray::benchmarks;
using namespace detray::algebra;

/// Run vector benchmarks
int main(int argc, char** argv) {
  constexpr std::size_t n_samples{100000};

  //
  // Prepare benchmarks
  //
  detray::benchmarks::configuration cfg_s{};
  // Reduce the number of samples, since a single SoA struct contains multiple
  // vectors
  cfg_s.n_samples(n_samples / Vc::float_v::Size);
  // Leave this to google benchmark
  cfg_s.do_warmup(false);

  // For double precision we need more samples (less vectors per SoA)
  detray::benchmarks::configuration cfg_d{cfg_s};
  cfg_d.n_samples(n_samples / Vc::double_v::Size);
  // Leave this to google benchmark
  cfg_d.do_warmup(false);

  using trf_f_t = transform3_bm<vc_soa::transform3<float>>;
  using trf_d_t = transform3_bm<vc_soa::transform3<double>>;

  using mat44_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                             bench_op::matrix::transpose>;
  using mat44_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                             bench_op::matrix::transpose>;
  using mat66_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                             bench_op::matrix::transpose>;
  using mat66_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                             bench_op::matrix::transpose>;
  using mat88_transp_f_t = matrix_unaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                             bench_op::matrix::transpose>;
  using mat88_transp_d_t = matrix_unaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                             bench_op::matrix::transpose>;

  using mat44_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                           bench_op::matrix::add>;
  using mat44_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                           bench_op::matrix::add>;
  using mat66_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                           bench_op::matrix::add>;
  using mat66_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                           bench_op::matrix::add>;
  using mat88_add_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                           bench_op::matrix::add>;
  using mat88_add_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                           bench_op::matrix::add>;

  using mat44_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 4, 4>,
                                           bench_op::matrix::mul>;
  using mat44_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 4, 4>,
                                           bench_op::matrix::mul>;
  using mat66_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 6, 6>,
                                           bench_op::matrix::mul>;
  using mat66_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 6, 6>,
                                           bench_op::matrix::mul>;
  using mat88_mul_f_t = matrix_binaryOP_bm<vc_soa::matrix_type<float, 8, 8>,
                                           bench_op::matrix::mul>;
  using mat88_mul_d_t = matrix_binaryOP_bm<vc_soa::matrix_type<double, 8, 8>,
                                           bench_op::matrix::mul>;

  using mat44_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 4, 4>,
                                         vc_soa::vector_type<float, 4>>;
  using mat44_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 4, 4>,
                                         vc_soa::vector_type<double, 4>>;
  using mat66_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 6, 6>,
                                         vc_soa::vector_type<float, 6>>;
  using mat66_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 6, 6>,
                                         vc_soa::vector_type<double, 6>>;
  using mat88_vec_f_t = matrix_vector_bm<vc_soa::matrix_type<float, 8, 8>,
                                         vc_soa::vector_type<float, 8>>;
  using mat88_vec_d_t = matrix_vector_bm<vc_soa::matrix_type<double, 8, 8>,
                                         vc_soa::vector_type<double, 8>>;

  std::cout << "-----------------------------------------------\n"
            << "Detray linear algebra benchmark (Vc SoA)\n"
            << "-----------------------------------------------\n\n"
            << "(single)\n"
            << cfg_s << "(double)\n"
            << cfg_d;

  //
  // Register all benchmarks
  //

  // Vector benchmarks
  DETRAY_DEFINE_VECTOR_BENCH(vc_soa)
  DETRAY_REGISTER_VECTOR_BENCH(cfg_s, cfg_d)

  // Transform benchmarks
  detray::benchmarks::register_benchmark<trf_f_t>(cfg_s, "_SINGLE");
  detray::benchmarks::register_benchmark<trf_d_t>(cfg_d, "_DOUBLE");

  // Matrix benchmarks
  detray::benchmarks::register_benchmark<mat44_transp_f_t>(cfg_s,
                                                           "_4x4_SINGLE");
  detray::benchmarks::register_benchmark<mat44_transp_d_t>(cfg_d,
                                                           "_4x4_DOUBLE");
  detray::benchmarks::register_benchmark<mat66_transp_f_t>(cfg_s,
                                                           "_6x6_SINGLE");
  detray::benchmarks::register_benchmark<mat66_transp_d_t>(cfg_d,
                                                           "_6x6_DOUBLE");
  detray::benchmarks::register_benchmark<mat88_transp_f_t>(cfg_s,
                                                           "_8x8_SINGLE");
  detray::benchmarks::register_benchmark<mat88_transp_d_t>(cfg_d,
                                                           "_8x8_DOUBLE");

  detray::benchmarks::register_benchmark<mat44_add_f_t>(cfg_s, "_4x4_SINGLE");
  detray::benchmarks::register_benchmark<mat44_add_d_t>(cfg_d, "_4x4_DOUBLE");
  detray::benchmarks::register_benchmark<mat66_add_f_t>(cfg_s, "_6x6_SINGLE");
  detray::benchmarks::register_benchmark<mat66_add_d_t>(cfg_d, "_6x6_DOUBLE");
  detray::benchmarks::register_benchmark<mat88_add_f_t>(cfg_s, "_8x8_SINGLE");
  detray::benchmarks::register_benchmark<mat88_add_d_t>(cfg_d, "_8x8_DOUBLE");

  detray::benchmarks::register_benchmark<mat44_mul_f_t>(cfg_s, "_4x4_SINGLE");
  detray::benchmarks::register_benchmark<mat44_mul_d_t>(cfg_d, "_4x4_DOUBLE");
  detray::benchmarks::register_benchmark<mat66_mul_f_t>(cfg_s, "_6x6_SINGLE");
  detray::benchmarks::register_benchmark<mat66_mul_d_t>(cfg_d, "_6x6_DOUBLE");
  detray::benchmarks::register_benchmark<mat88_mul_f_t>(cfg_s, "_8x8_SINGLE");
  detray::benchmarks::register_benchmark<mat88_mul_d_t>(cfg_d, "_8x8_DOUBLE");

  detray::benchmarks::register_benchmark<mat44_vec_f_t>(cfg_s, "_4x4_SINGLE");
  detray::benchmarks::register_benchmark<mat44_vec_d_t>(cfg_d, "_4x4_DOUBLE");
  detray::benchmarks::register_benchmark<mat66_vec_f_t>(cfg_s, "_6x6_SINGLE");
  detray::benchmarks::register_benchmark<mat66_vec_d_t>(cfg_d, "_6x6_DOUBLE");
  detray::benchmarks::register_benchmark<mat88_vec_f_t>(cfg_s, "_8x8_SINGLE");
  detray::benchmarks::register_benchmark<mat88_vec_d_t>(cfg_d, "_8x8_DOUBLE");

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
