// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).
#include <detray/definitions/units.hpp>

namespace traccc {

template <typename scalar_t>
using unit = detray::unit<scalar_t>;

template <typename scalar_t>
using constant = detray::constant<scalar_t>;

// epsilon for float variables
constexpr scalar float_epsilon = 1e-5f;

}  // namespace traccc
