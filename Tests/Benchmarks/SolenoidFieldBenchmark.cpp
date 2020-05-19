// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <string>

using namespace Acts::UnitLiterals;

int main(int argc, char* argv[]) {
  size_t iters_map = 5e2;
  size_t iters_solenoid = 3;
  size_t runs_solenoid = 1000;
  if (argc >= 2) {
    iters_map = std::stoi(argv[1]);
  }
  if (argc >= 3) {
    iters_solenoid = std::stoi(argv[2]);
  }
  if (argc >= 4) {
    runs_solenoid = std::stoi(argv[3]);
  }

  const double L = 5.8_m;
  const double R = (2.56 + 2.46) * 0.5 * 0.5_m;
  const size_t nCoils = 1154;
  const double bMagCenter = 2_T;
  const size_t nBinsR = 150;
  const size_t nBinsZ = 200;

  double rMin = -0.1;
  double rMax = R * 2.;
  double zMin = 2 * (-L / 2.);
  double zMax = 2 * (L / 2.);

  Acts::SolenoidBField bSolenoidField({R, L, nCoils, bMagCenter});
  std::cout << "Building interpolated field map" << std::endl;
  auto mapper = Acts::solenoidFieldMapper({rMin, rMax}, {zMin, zMax},
                                          {nBinsR, nBinsZ}, bSolenoidField);
  using BField_t = Acts::InterpolatedBFieldMap<decltype(mapper)>;

  BField_t::Config cfg(std::move(mapper));
  auto bFieldMap = BField_t(std::move(cfg));
  using Cache = typename BField_t::Cache;
  Acts::MagneticFieldContext mctx{};

  std::minstd_rand rng;
  std::uniform_real_distribution<> zDist(1.5 * (-L / 2.), 1.5 * L / 2.);
  std::uniform_real_distribution<> rDist(0, R * 1.5);
  std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
  auto genPos = [&]() -> Acts::Vector3D {
    const double z = zDist(rng), r = rDist(rng), phi = phiDist(rng);
    return {r * std::cos(phi), r * std::sin(phi), z};
  };

  auto csv = [&](const std::string& name, auto res) {

    auto& os = std::cerr;

    os << name << "," << res.run_timings.size() << "," 
      << res.iters_per_run << ","
      << res.totalTime().count() << ","
      << res.runTimeMedian().count() << ","
      << 1.96*res.runTimeError().count() << ","
      << res.iterTimeAverage().count() << ","
      << 1.96*res.iterTimeError().count();

    os << std::endl;
  };

  std::cerr << "name,runs,iters,total_time,run_time_median,run_time_error,iter_time_average,iter_time_error" << std::endl;

  // SolenoidBField lookup is so slow that the cost of generating a random field
  // lookup position is negligible in comparison...
  std::cout << "Benchmarking random SolenoidBField lookup: " << std::flush;
  const auto solenoid_result = Acts::Test::microBenchmark(
      [&] { return bSolenoidField.getField(genPos()); }, iters_solenoid,
      runs_solenoid);
  std::cout << solenoid_result << std::endl;
  csv("solenoid", solenoid_result);

  // ...but for interpolated B-field map, the overhead of a field lookup is
  // comparable to that of generating a random position, so we must be more
  // careful. Hence we do two microbenchmarks which represent a kind of
  // lower and upper bound on field lookup performance.
  //
  // - The first benchmark operates at constant position, so it measures only
  //   field lookup overhead but has unrealistically good cache locality. In
  //   that sense, it provides a lower bound of field lookup performance.
  std::cout << "Benchmarking interpolated field lookup: " << std::flush;
  const auto fixedPos = genPos();
  const auto map_cached_result = Acts::Test::microBenchmark(
      [&] { return bFieldMap.getField(fixedPos); }, iters_map);
  std::cout << map_cached_result << std::endl;
  csv("interp_nocache_fixed", map_cached_result);


  // - The second benchmark generates random positions, so it is biased by the
  //   cost of random position generation and has unrealistically bad cache
  //   locality, but provides an upper bound of field lookup performance.
  std::cout << "Benchmarking random interpolated field lookup: " << std::flush;
  const auto map_rand_result = Acts::Test::microBenchmark(
      [&] { return bFieldMap.getField(genPos()); }, iters_map);
  std::cout << map_rand_result << std::endl;
  csv("interp_nocache_random", map_rand_result);


  {
  std::cout << "Benchmarking cached interpolated field lookup: " << std::flush;
  Cache cache{mctx};
  const auto map_cached_result_cache = Acts::Test::microBenchmark(
      [&] { return bFieldMap.getField(fixedPos, cache); }, iters_map);
  std::cout << map_cached_result_cache << std::endl;
  csv("interp_cache_fixed", map_cached_result_cache);
}


  // - The second benchmark generates random positions, so it is biased by the
  //   cost of random position generation and has unrealistically bad cache
  //   locality, but provides an upper bound of field lookup performance.
  {
    std::cout << "Benchmarking cached random interpolated field lookup: " << std::flush;
    Cache cache2{mctx};
    const auto map_rand_result_cache = Acts::Test::microBenchmark(
        [&] { return bFieldMap.getField(genPos(), cache2); }, iters_map);
    std::cout << map_rand_result_cache << std::endl;
    csv("interp_cache_random", map_rand_result_cache);
  }

  {
    std::cout << "Benchmarking advancing interpolated field lookup: " << std::flush;
    Acts::Vector3D pos{0,0,0};
    Acts::Vector3D dir{};
    dir.setRandom();
    double h = 1e-3;
    const auto map_adv_result = Acts::Test::microBenchmark(
        [&] {
        pos += dir * h;
        return bFieldMap.getField(pos); 
        }, iters_map);
    std::cout << map_adv_result << std::endl;
    csv("interp_nocache_adv", map_adv_result);
  }  

  {
    std::cout << "Benchmarking cached advancing interpolated field lookup: " << std::flush;
    Cache cache{mctx};
    Acts::Vector3D pos{0,0,0};
    Acts::Vector3D dir{};
    dir.setRandom();
    double h = 1e-3;
    const auto map_adv_result_cache = Acts::Test::microBenchmark(
        [&] {
        pos += dir * h;
        return bFieldMap.getField(pos, cache); 
        }, iters_map);
    std::cout << map_adv_result_cache << std::endl;
    csv("interp_cache_adv", map_adv_result_cache);
  }  
}
