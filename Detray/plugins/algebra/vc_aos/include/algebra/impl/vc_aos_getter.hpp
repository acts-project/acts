// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/matrix_getter.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::vc_aos::storage {

using algebra::storage::block;
using algebra::storage::element;
using algebra::storage::set_block;

/// Get a vector of a const matrix
template <std::size_t SIZE, std::size_t ROW, std::size_t COL,
          concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr decltype(auto) vector(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL> &m,
    const std::size_t row, const std::size_t col) noexcept {
  return algebra::storage::block_getter{}.template vector<SIZE>(m, row, col);
}

}  // namespace detray::algebra::vc_aos::storage
