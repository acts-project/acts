// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/propagator/codegen/transport_jacobian.hpp"

namespace detray {

namespace concepts {

template <typename T>
concept transport_jacobian =
    requires { typename T::algebra_type; } &&
    (std::same_as<T, detray::detail::transport_jacobian_matrix_with_gradient<
                         typename T::algebra_type>> ||
     std::same_as<T, detray::detail::transport_jacobian_matrix_without_gradient<
                         typename T::algebra_type>>);

}  // namespace concepts

namespace matrix {

template <concepts::transport_jacobian matrix_type>
DETRAY_HOST_DEVICE constexpr auto identity() {
  return matrix_type::identity();
}

}  // namespace matrix

namespace traits {

template <typename algebra_t>
struct dimensions<
    ::detray::detail::transport_jacobian_matrix_with_gradient<algebra_t>> {
  using size_type = std::size_t;

  static constexpr size_type _dim{2};
  static constexpr size_type _rows{8};
  static constexpr size_type _columns{8};
};

template <typename algebra_t>
struct dimensions<
    ::detray::detail::transport_jacobian_matrix_without_gradient<algebra_t>> {
  using size_type = std::size_t;

  static constexpr size_type _dim{2};
  static constexpr size_type _rows{8};
  static constexpr size_type _columns{8};
};

template <typename algebra_t>
struct index<
    ::detray::detail::transport_jacobian_matrix_with_gradient<algebra_t>> {
  using type = std::size_t;
};

template <typename algebra_t>
struct index<
    ::detray::detail::transport_jacobian_matrix_without_gradient<algebra_t>> {
  using type = std::size_t;
};

}  // namespace traits
}  // namespace detray
