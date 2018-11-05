// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/BFieldMapUtils.hpp"
#include "Acts/Utilities/Units.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <string>

int
main(int argc, char* argv[])
{

  size_t n = 1e6;
  if (argc == 2) {
    n = std::stoi(argv[1]);
  }

  const double L          = 5.8 * Acts::units::_m;
  const double R          = (2.56 + 2.46) * 0.5 * 0.5 * Acts::units::_m;
  const size_t nCoils     = 1154;
  const double bMagCenter = 2. * Acts::units::_T;
  const size_t nBinsR     = 150;
  const size_t nBinsZ     = 200;

  double rMin = -0.1;
  double rMax = R * 2.;
  double zMin = 2 * (-L / 2.);
  double zMax = 2 * (L / 2.);

  size_t printStep = n / 10;

  Acts::SolenoidBField bSolenoidField({R, L, nCoils, bMagCenter});
  std::cout << "building map" << std::endl;
  auto mapper = Acts::solenoidFieldMapper(
      {rMin, rMax}, {zMin, zMax}, {nBinsR, nBinsZ}, bSolenoidField);
  Acts::InterpolatedBFieldMap::Config cfg;
  cfg.mapper     = std::move(mapper);
  auto bFieldMap = Acts::InterpolatedBFieldMap(std::move(cfg));

  std::mt19937                     rng;
  std::uniform_real_distribution<> zDist(1.5 * (-L / 2.), 1.5 * L / 2.);
  std::uniform_real_distribution<> rDist(0, R * 1.5);
  std::uniform_real_distribution<> phiDist(-M_PI, M_PI);

  std::cout << "number of points: " << n << std::endl;
  std::cout << "start" << std::endl;
  using clock = std::chrono::steady_clock;
  auto start  = clock::now();

  double         z, r, phi;
  Acts::Vector3D pos;
  Acts::Vector3D B;
  for (size_t i = 0; i < n; i++) {
    if (i % printStep == 0) {
      std::cout << i << std::endl;
    }

    z   = zDist(rng);
    r   = rDist(rng);
    phi = phiDist(rng);
    pos = {r * std::cos(phi), r * std::sin(phi), z};

    B = bSolenoidField.getField(pos);
  }

  auto   end = clock::now();
  double ms_solenoid
      = std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
            .count();
  start = clock::now();

  for (size_t i = 0; i < n; i++) {
    if (i % printStep == 0) {
      std::cout << i << std::endl;
    }

    z   = zDist(rng);
    r   = rDist(rng);
    phi = phiDist(rng);
    pos = {r * std::cos(phi), r * std::sin(phi), z};

    B = bFieldMap.getField(pos);
  }

  end = clock::now();
  auto ms_map
      = std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
            .count();

  std::cout << "solenoid: " << (ms_solenoid * 1000.)
            << "ms, per lookup: " << (ms_solenoid / n * 1000.) << "ms"
            << std::endl;
  std::cout << "map: " << (ms_map * 1000.)
            << "ms, per lookup: " << (ms_map / n * 1000.) << "ms" << std::endl;
}
