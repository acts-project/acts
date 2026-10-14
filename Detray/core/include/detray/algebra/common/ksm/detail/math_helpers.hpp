// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"
#include "detray/algebra/common/ksm/detail/symbolic.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <concepts>

// This file contains a few helper functions, in particular related to
// matrix determinants.
namespace detray::ksm::detail {
/// @brief Helper functions for getting a matrix determinant as a symbolic
/// value.
///
/// The only important thing this function can return is whether the
/// determinant is structurally 0, structurally non-zero, or variable.
///
/// @{
template <typename T>
struct get_determinant_helper {};

// Definitionally, the determinant of a 1x1 matrix is the only value in that
// matrix.
template <typename T>
  requires(T::rows == 1 && T::columns == 1)
struct get_determinant_helper<T> {
  using type = T::template value_at<0, 0>;
};

// As we know, the determinant of [[a, b], [c, d]] is ad-bc. We compute that
// symbolically here and resolve it to a value.
template <typename T>
  requires(T::rows == 2 && T::columns == 2)
struct get_determinant_helper<T> {
 private:
  using t1 = value_multiplication<typename T::template value_at<0, 0>,
                                  typename T::template value_at<1, 1>>;
  using t2 = value_multiplication<typename T::template value_at<1, 0>,
                                  typename T::template value_at<0, 1>>;

  // Recall the syntax of resolve_value: the first bool asks whether all
  // values are from the same matrix; the second asks if we should collapse
  // the result to a variable. We don't want to do that yet, just in case
  // ad and bc are structurally identical. If they are, collapsing to
  // variables would yield V1-V2 which is not structurally zero.
  using r1 = resolve_value<t1, true, false>::type;
  using r2 = resolve_value<t2, true, false>::type;

 public:
  // Then model the symbolic subtraction. This we _do_ want to resolve to a
  // variable.
  using type = resolve_value<value_subtraction<r1, r2>, true, true>::type;
};
/// @}
}  // namespace detray::ksm::detail

namespace detray::ksm::concepts {
/// @brief True iff the matrix is possibly invertible, i.e. if the determinant
/// is not structurally zero.
///
/// @warning This doesn't guarantee that a matrix is invertible at runtime, as
/// the runtime values might still make inversion impossible.
template <typename T>
concept is_invertible =
    is_canonical<T> && T::rows == T::columns &&
    (T::rows == 1 || T::rows == 2) &&
    !std::same_as<typename detail::get_determinant_helper<T>::type, zero>;
}  // namespace detray::ksm::concepts
