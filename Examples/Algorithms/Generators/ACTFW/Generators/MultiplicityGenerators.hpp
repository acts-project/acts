// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <random>

#include "ACTFW/Framework/RandomNumbers.hpp"

namespace FW {

struct FixedMultiplicityGenerator {
  size_t n = 1;

  size_t operator()(RandomEngine& /* unused */) const { return n; }
};

struct PoissonMultiplicityGenerator {
  size_t mu = 1;

  size_t operator()(RandomEngine& rng) const {
    return (0 < mu) ? std::poisson_distribution<size_t>(mu)(rng) : 0;
  }
};

}  // namespace FW
