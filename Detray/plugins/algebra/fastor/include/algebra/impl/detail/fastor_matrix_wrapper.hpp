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

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

// System include(s).
#include <cstddef>

namespace detray::algebra::fastor {

/// @brief Fastor Matrix type
///
/// This class is needed because Fastor differentiates between normal
/// matrix-matrix multiplication and element-wise matrix-matrix multiplication.
///
/// In Fastor, the former is performed by using `operator%` or the
/// `Fastor::matmul()` function, whereas the latter is done using `operator*`.
/// However, the algebra-plugins repository expects operator* for normal
/// matrix-matrix multiplication. To resolve this issue, this wrapper class
/// around `Fastor::Tensor` was created.
///
/// The class inherits from `Fastor::Tensor` because we want objects of this
/// class to behave the way a `Fastor::Tensor` would. Inheriting from
/// `Fastor::Tensor` allows this class to reuse all the functions defined in the
/// parent class (i.e. `Fastor::Tensor`).
template <concepts::scalar T, std::size_t M1, std::size_t N>
class Matrix : public Fastor::Tensor<T, M1, N> {
 public:
  /// Inherit all constructors from the base class
  using Fastor::Tensor<T, M1, N>::Tensor;

 private:
  /// When we encounter an `operator*` function call between Matrix objects,
  /// we will catch it and handle it correctly by invoking `Fastor::matmul()`.
  ///
  /// The data type contained in the `other` has a separate template parameter
  /// dedicated to it because in certain cases, we might want to multiply say,
  /// a float matrix with a double matrix and not have it produce a
  /// compilation error.
  ///
  /// The `static_cast` is there to signal both to the compiler and the reader
  /// that we wish to interpret the `Matrix` object as a `Fastor::Tensor`
  /// here.
  /// @{
  template <std::size_t LR, std::size_t C, std::size_t RC, concepts::scalar S>
  constexpr friend Matrix<S, LR, RC> operator*(const Matrix<S, LR, C>& lhs,
                                               const Matrix<S, C, RC>& rhs);

  template <std::size_t R, std::size_t C, concepts::scalar S>
  constexpr friend bool operator==(const Matrix<S, R, C>& lhs,
                                   const Matrix<S, R, C>& rhs);
  /// @}

};  // class Matrix

template <std::size_t LR, std::size_t C, std::size_t RC, concepts::scalar S>
constexpr Matrix<S, LR, RC> operator*(const Matrix<S, LR, C>& lhs,
                                      const Matrix<S, C, RC>& rhs) {
  return Fastor::matmul(static_cast<Fastor::Tensor<S, LR, C>>(lhs),
                        static_cast<Fastor::Tensor<S, C, RC>>(rhs));
}

template <std::size_t R, std::size_t C, concepts::scalar S>
constexpr bool operator==(const Matrix<S, R, C>& lhs,
                          const Matrix<S, R, C>& rhs) {
  return Fastor::isequal(static_cast<Fastor::Tensor<S, R, C>>(lhs),
                         static_cast<Fastor::Tensor<S, R, C>>(rhs), 0.f);
}

}  // namespace detray::algebra::fastor
