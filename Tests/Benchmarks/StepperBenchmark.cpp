// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

#include <memory>

#include "StepperBenchmarkCommons.hpp"

using namespace Acts;
using namespace Acts::Test;

using Stepper = EigenStepper<>;

int main(int argc, char* argv[]) {
  BenchmarkStepper benchmark;
  if (auto ret = benchmark.parseOptions(argc, argv)) {
    return *ret;
  }
  std::shared_ptr<MagneticFieldProvider> bField = benchmark.makeField();
  AtlasStepper atlasStepper(bField);
  benchmark.run(atlasStepper, "AtlasStepper");
  EigenStepper eigenStepper(bField);
  benchmark.run(eigenStepper, "EigenStepper");
  StraightLineStepper straightLineStepper;
  benchmark.run(straightLineStepper, "StraightLineStepper");
  return 0;
}
