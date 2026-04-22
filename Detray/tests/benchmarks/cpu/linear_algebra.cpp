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
#include "detray/benchmarks/types.hpp"
#include "detray/benchmarks/cpu/matrix_benchmark.hpp"
#include "detray/benchmarks/cpu/vector_benchmark.hpp"
#include "detray/benchmarks/cpu/transform_benchmark.hpp"
// clang-format on

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>

using namespace detray::benchmarks;

/// Run vector benchmarks
int main(int argc, char** argv) {
  //
  // Prepare benchmarks
  //
  benchmark_base::configuration cfg{};
  cfg.n_samples(100000);
  // Leave this to google benchmark
  cfg.do_warmup(false);

  std::cout << "-----------------------------------------------\n"
            << "Detray linear algebra benchmark ("
            << detray::benchmarks::plugin_name << ")\n"
            << "-----------------------------------------------\n\n"
            << cfg;

//
// Define and register all benchmarks
//
#if DETRAY_ALGEBRA_ARRAY
  DETRAY_DEFINE_VECTOR_BENCH(detray::algebra::array)
  DETRAY_DEFINE_TRANSFORM_BENCH(detray::algebra::array)
  DETRAY_DEFINE_MATRIX_BENCH(detray::algebra::array)
#elif DETRAY_ALGEBRA_EIGEN
  DETRAY_DEFINE_VECTOR_BENCH(detray::algebra::eigen)
  DETRAY_DEFINE_TRANSFORM_BENCH(detray::algebra::eigen)
  DETRAY_DEFINE_MATRIX_BENCH(detray::algebra::eigen)
#elif DETRAY_ALGEBRA_FASTOR
  DETRAY_DEFINE_VECTOR_BENCH(detray::algebra::fastor)
  DETRAY_DEFINE_TRANSFORM_BENCH(detray::algebra::fastor)
  DETRAY_DEFINE_MATRIX_BENCH(detray::algebra::fastor)
#elif DETRAY_ALGEBRA_SMATRIX
  DETRAY_DEFINE_VECTOR_BENCH(detray::algebra::smatrix)
  DETRAY_DEFINE_TRANSFORM_BENCH(detray::algebra::smatrix)
  DETRAY_DEFINE_MATRIX_BENCH(detray::algebra::smatrix)
#elif DETRAY_ALGEBRA_VC_AOS
  DETRAY_DEFINE_VECTOR_BENCH(detray::algebra::vc_aos)
  DETRAY_DEFINE_TRANSFORM_BENCH(detray::algebra::vc_aos)
  DETRAY_DEFINE_MATRIX_BENCH(detray::algebra::vc_aos)
#endif

  DETRAY_REGISTER_VECTOR_BENCH(cfg, cfg)
  DETRAY_REGISTER_TRANSFORM_BENCH(cfg, cfg)
  DETRAY_REGISTER_MATRIX_BENCH(cfg, cfg)

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
