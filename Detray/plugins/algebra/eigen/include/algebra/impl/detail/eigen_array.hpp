// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20012
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
// Suppress GCC false positives in Eigen while keeping our own code checked.
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && \
    __GNUC__ >= 12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Core>
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && \
    __GNUC__ >= 12
#pragma GCC diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

// System include(s).
#include <cstddef>

namespace detray::algebra::eigen {

/// Eigen array type
template <concepts::scalar T, int N>
class array : public Eigen::Matrix<T, N, 1, 0, N, 1> {
 public:
  /// Inherit all constructors from the base class
  using Eigen::Matrix<T, N, 1, 0, N, 1>::Matrix;

};  // class array

}  // namespace detray::algebra::eigen
