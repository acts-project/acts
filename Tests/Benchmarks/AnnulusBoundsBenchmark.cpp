// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

using namespace Acts;

int main(int /*argc*/, char** /*argv[]*/) {
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> xDist(-2, 20);
  std::uniform_real_distribution<double> yDist(-2, 20);

  std::ofstream os{"annulus.csv"};
  os << "x,y,inside_abs,inside_tol0,inside_tol1,inside_tol01,inside_cov"
     << std::endl;

  // === PROBLEM DATA ===

  // bounds definition
  double minRadius = 7.2;
  double maxRadius = 12.0;
  double minPhi = 0.74195;
  double maxPhi = 1.33970;
  Vector2 offset(-2., 2.);
  AnnulusBounds aBounds(minRadius, maxRadius, minPhi, maxPhi, offset);

  // helper to convert to expected local frame
  auto toStripFrame = [&](const Vector2& xy) -> Vector2 {
    auto shifted = xy + offset;
    double r = VectorHelpers::perp(shifted);
    double phi = VectorHelpers::phi(shifted);
    return Vector2(r, phi);
  };

  auto random_point = [&]() -> Vector2 { return {xDist(rng), yDist(rng)}; };

  // for covariance based check, set up one;
  ActsMatrix<2, 2> cov;
  cov << 1.0, 0, 0, 0.05;

  BoundaryCheck bcAbs{true};
  BoundaryCheck bcTol0{true, false, 1.0, 0};
  BoundaryCheck bcTol1{false, true, 0, 0.2};
  BoundaryCheck bcTol01{true, true, 1.0, 0.2};
  BoundaryCheck bcCov{cov, 1};

  // visualization to make sense of things
  for (size_t i = 0; i < 10000; i++) {
    const Vector2 loc{xDist(rng), yDist(rng)};
    auto locPC = toStripFrame(loc);
    bool isInsideAbs = aBounds.inside(locPC, bcAbs);
    bool isInsideTol0 = aBounds.inside(locPC, bcTol0);
    bool isInsideTol1 = aBounds.inside(locPC, bcTol1);
    bool isInsideTol01 = aBounds.inside(locPC, bcTol01);

    bool isInsideCov = aBounds.inside(locPC, bcCov);

    os << loc.x() << "," << loc.y() << "," << isInsideAbs << "," << isInsideTol0
       << "," << isInsideTol1 << "," << isInsideTol01 << "," << isInsideCov
       << std::endl;
  }

  std::vector<std::tuple<Vector2, bool, std::string>> testPoints{{
      {{7.0, 6.0}, true, "center"},
      {{3.0, 1.0}, false, "radial far out low"},
      {{5.0, 4.0}, false, "radial close out low"},
      {{6.0, 5.0}, true, "radial close in low"},
      {{8.5, 7.5}, true, "radial close in high"},
      {{10.0, 9.0}, false, "radial close out high"},
      {{15.0, 15.0}, false, "radial far out high"},
      {{0.0, 11.0}, false, "angular far out left"},
      {{4.0, 9.0}, false, "angular close out left"},
      {{5.0, 8.5}, true, "angular close in left"},
      {{9.0, 5.0}, true, "angular close in right"},
      {{9.5, 4.5}, false, "angular close out right"},
      {{11.0, 0.0}, false, "angular far out right"},
      {{2.0, 4.0}, false, "quad low left"},
      {{5.0, 14.0}, false, "quad high left"},
      {{13.0, 6.0}, false, "quad high right"},
      {{7.0, 1.0}, false, "quad low right"},
  }};

  // === BENCHMARKS ===

  // Number of benchmark runs
  constexpr int NTESTS = 5'000;

  // Some checks are much slower, so we tune down benchmark iterations
  constexpr int NTESTS_SLOW = NTESTS / 10;

  // We use this to switch between iteration counts
  enum class Mode { FastOutside, SlowOutside };

  // Benchmark output display
  auto print_bench_header = [](const std::string& check_name) {
    std::cout << check_name << ":" << std::endl;
  };
  auto print_bench_result = [](const std::string& bench_name,
                               const Acts::Test::MicroBenchmarkResult& res) {
    std::cout << "- " << bench_name << ": " << res << std::endl;
  };

  // Benchmark runner
  auto run_bench = [&](auto&& iteration, int num_iters,
                       const std::string& bench_name) {
    auto bench_result = Acts::Test::microBenchmark(iteration, num_iters);
    print_bench_result(bench_name, bench_result);
  };

  auto run_bench_with_inputs = [&](auto&& iterationWithArg, auto&& inputs,
                                   const std::string& bench_name) {
    auto bench_result = Acts::Test::microBenchmark(iterationWithArg, inputs);
    print_bench_result(bench_name, bench_result);
  };
  auto run_all_benches = [&](const BoundaryCheck& check,
                             const std::string& check_name, const Mode mode) {
    // Announce a set of benchmarks
    print_bench_header(check_name);

    // Pre-determined "interesting" test points
    int num_inside_points = 0;
    int num_outside_points = 0;
    switch (mode) {
      case Mode::FastOutside:
        num_inside_points = NTESTS;
        num_outside_points = NTESTS;
        break;
      case Mode::SlowOutside:
        num_inside_points = NTESTS;
        num_outside_points = NTESTS_SLOW;
    };

    for (const auto& [loc, inside, label] : testPoints) {
      const auto locPC = toStripFrame(loc);
      run_bench([&] { return aBounds.inside(locPC, check); },
                inside ? num_inside_points : num_outside_points, label);
    }

    // Pre-rolled random points
    std::vector<Vector2> points(num_outside_points);
    std::generate(points.begin(), points.end(),
                  [&] { return toStripFrame(random_point()); });
    run_bench_with_inputs(
        [&](const auto& point) { return aBounds.inside(point, check); }, points,
        "Random");
  };

  // Benchmark scenarios
  run_all_benches(bcAbs, "Absolute", Mode::FastOutside);
  run_all_benches(bcTol0, "Tolerance 0", Mode::FastOutside);
  run_all_benches(bcTol1, "Tolerance 1", Mode::FastOutside);
  run_all_benches(bcTol01, "Tolerance 01", Mode::FastOutside);
  run_all_benches(bcCov, "Covariance", Mode::SlowOutside);

  return 0;
}
