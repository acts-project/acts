// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

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

/// @brief Fastor vector type that provides equality op
template <concepts::scalar T, std::size_t N>
class Vector : public Fastor::Tensor<T, N> {
 public:
  /// Inherit all constructors from the base class
  using Fastor::Tensor<T, N>::Tensor;

 private:
  /// Multiplication
  template <std::size_t C, concepts::scalar S, typename D>
  constexpr friend auto operator*(const Fastor::AbstractTensor<D, 2>& lhs,
                                  const Vector<S, C>& vector);

  template <std::size_t C, std::size_t R, concepts::scalar S>
  constexpr friend Vector<S, R> operator*(const Fastor::Tensor<S, R, C>& lhs,
                                          const Vector<S, C>& vector);

  /// Equality operator for fastor vectors
  template <std::size_t M, concepts::scalar S>
  constexpr friend bool operator==(const Vector<S, M>& lhs,
                                   const Vector<S, M>& rhs);
  /// @}

};  // class Vector

template <std::size_t C, concepts::scalar S, typename D>
constexpr auto operator*(const Fastor::AbstractTensor<D, 2>& lhs,
                         const Vector<S, C>& vector) {
  return Fastor::matmul(Fastor::evaluate(lhs),
                        static_cast<Fastor::Tensor<S, C>>(vector));
}

template <std::size_t C, std::size_t R, concepts::scalar S>
constexpr Vector<S, R> operator*(const Fastor::Tensor<S, R, C>& lhs,
                                 const Vector<S, C>& vector) {
  return Fastor::matmul(lhs, static_cast<Fastor::Tensor<S, C>>(vector));
}

template <std::size_t M, concepts::scalar S>
constexpr bool operator==(const Vector<S, M>& lhs, const Vector<S, M>& rhs) {
  return Fastor::isequal(static_cast<Fastor::Tensor<S, M>>(lhs),
                         static_cast<Fastor::Tensor<S, M>>(rhs), 0.f);
}

}  // namespace detray::algebra::fastor
