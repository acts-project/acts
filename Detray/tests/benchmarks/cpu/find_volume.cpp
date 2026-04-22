// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/utils/logging.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Google include(s).
#include <benchmark/benchmark.h>

// System include(s).
#include <iostream>

// Use the detray:: namespace implicitly.
using namespace detray;

using bench_algebra = benchmarks::algebra;
using scalar = benchmarks::scalar;

// Benchmarks the cost of searching a volume by position
void BM_FIND_VOLUMES(benchmark::State &state) {
  // This test is broken at the moment, don't run it.
  state.SkipWithError("Benchmark disabled for the toy geometry");
  return;

  // Detector configuration
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.n_edc_layers(7u);
  auto [d, names] = build_toy_detector<bench_algebra>(host_mr, toy_cfg);

  static const unsigned int itest = 10000u;

  constexpr auto vol_grid_idx{decltype(d)::accel::id::e_volume_default};
  auto volume_grid = d.accelerator_store().template get<vol_grid_idx>()[0];

  const auto &axis_r = volume_grid.get_axis<axis::label::e_r>();
  const auto &axis_z = volume_grid.get_axis<axis::label::e_z>();

  // Get a rough step size from irregular axes
  auto range0 = axis_r.span();
  auto range1 = axis_z.span();

  scalar step0{(range0[1] - range0[0]) / itest};
  scalar step1{(range1[1] - range1[0]) / itest};

  std::size_t successful{0u};
  std::size_t unsuccessful{0u};

  for (auto _ : state) {
    for (unsigned int i1 = 0u; i1 < itest; ++i1) {
      for (unsigned int i0 = 0u; i0 < itest; ++i0) {
        benchmarks::vector3 rz{static_cast<scalar>(i0) * step0,
                               static_cast<scalar>(0.f),
                               static_cast<scalar>(i1) * step1};
        const auto &v = d.volume(rz);

        benchmark::DoNotOptimize(successful);
        benchmark::DoNotOptimize(unsuccessful);
        if (v.index() == dindex_invalid) {
          ++unsuccessful;
        } else {
          ++successful;
        }
        benchmark::ClobberMemory();
      }
    }
  }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST("Successful   : " << successful);
  DETRAY_INFO_HOST("Unsuccessful : " << unsuccessful);
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_FIND_VOLUMES)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
