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

  std::minstd_rand rng;
  std::uniform_real_distribution<> zDist(1.5 * (-L / 2.), 1.5 * L / 2.);
  std::uniform_real_distribution<> rDist(0, R * 1.5);
  std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
  auto genPos = [&]() -> Acts::Vector3D {
    const double z = zDist(rng), r = rDist(rng), phi = phiDist(rng);
    return {r * std::cos(phi), r * std::sin(phi), z};
  };

  // SolenoidBField lookup is so slow that the cost of generating a random field
  // lookup position is negligible in comparison...
  std::cout << "Benchmarking random SolenoidBField lookup: " << std::flush;
  const auto solenoid_result = Acts::Test::microBenchmark(
      [&] { return bSolenoidField.getField(genPos()); }, iters_solenoid,
      runs_solenoid);
  std::cout << solenoid_result << std::endl;

  // ...but for interpolated B-field map, the overhead of a field lookup is
  // comparable to that of generating a random position, so we must be more
  // careful. Hence we do two microbenchmarks which represent a kind of
  // lower and upper bound on field lookup performance.
  //
  // - The first benchmark operates at constant position, so it measures only
  //   field lookup overhead but has unrealistically good cache locality. In
  //   that sense, it provides a lower bound of field lookup performance.
  std::cout << "Benchmarking cached interpolated field lookup: " << std::flush;
  const auto fixedPos = genPos();
  const auto map_cached_result = Acts::Test::microBenchmark(
      [&] { return bFieldMap.getField(fixedPos); }, iters_map);
  std::cout << map_cached_result << std::endl;

  // - The second benchmark generates random positions, so it is biased by the
  //   cost of random position generation and has unrealistically bad cache
  //   locality, but provides an upper bound of field lookup performance.
  std::cout << "Benchmarking random interpolated field lookup: " << std::flush;
  const auto map_rand_result = Acts::Test::microBenchmark(
      [&] { return bFieldMap.getField(genPos()); }, iters_map);
  std::cout << map_rand_result << std::endl;
}
