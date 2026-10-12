// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/track_parametrization.hpp"

namespace traccc {

struct track_params_estimation_config {
  std::array<scalar, e_bound_size> initial_sigma = {
      1.f * unit<scalar>::mm,
      1.f * unit<scalar>::mm,
      1.f * unit<scalar>::degree,
      1.f * unit<scalar>::degree,
      0.f * unit<scalar>::e / unit<scalar>::GeV,
      1.f * unit<scalar>::ns};

  scalar initial_sigma_qopt = 0.1f * unit<scalar>::e / unit<scalar>::GeV;

  scalar initial_sigma_pt_rel = 0.1f;

  std::array<scalar, e_bound_size> initial_inflation = {1.f, 1.f, 1.f,
                                                        1.f, 1.f, 100.f};
};

}  // namespace traccc
