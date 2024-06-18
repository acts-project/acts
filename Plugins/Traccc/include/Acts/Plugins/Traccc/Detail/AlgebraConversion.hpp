// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"

// System include(s)
#include <cstddef>

namespace Acts::TracccPlugin::detail {

/// @brief Creates a new Acts vector from another vector type.
template <std::size_t N, typename dvector_t>
inline Acts::ActsVector<N> toActsVector(const dvector_t& dvec) {
  Acts::ActsVector<N> res = Eigen::Map<const Eigen::Matrix<typename dvector_t::value_type, N, 1>>(dvec.data()).template cast<ActsScalar>();
  return res;
}

/// @brief Creates a new Acts square matrix from another square matrix type.
template <std::size_t N, typename matrixNxN_t>
inline Acts::ActsSquareMatrix<N> toActsSquareMatrix(const matrixNxN_t& mat) {
  Acts::ActsSquareMatrix<N> res = Eigen::Map<const Eigen::Matrix<typename matrixNxN_t::value_type, N, N>>(mat.data()).template cast<ActsScalar>();
  return res;
}

}  // namespace Acts::TracccPlugin::detail
