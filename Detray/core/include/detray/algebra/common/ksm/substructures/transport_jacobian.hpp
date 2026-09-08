// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/known_substructure_matrix.hpp"

// System include(s)
#include <type_traits>

namespace detray::ksm {
/// @brief Substructure of the RK transport Jacobian.
///
/// The 8x8 free-to-free Jacobian the Runge-Kutta stepper accumulates. Both
/// the transport Jacobian and the update matrix have the same substructure,
/// so it can be accumulated (as long as it is closed under multiplication).
///
/// If the field gradient is enabled, the substructure is:
///
///      px py pz  t  dx dy dz qop
///  px [ v  v  v  0 | v  v  v | v ]
///  py [ v  v  v  0 | v  v  v | v ]
///  pz [ v  v  v  0 | v  v  v | v ]
///   t [ 0  0  0  1 | 0  0  0 | 0 ]
///  dx [ v  v  v  0 | v  v  v | v ]
///  dy [ v  v  v  0 | v  v  v | v ]
///  dz [ v  v  v  0 | v  v  v | v ]
/// qop [ 0  0  0  0 | 0  0  0 | v ]
///
/// Otherwise, it is:
///
///      px py pz  t  dx dy dz qop
///  px [ 1  0  0  0 | v  v  v | v ]
///  py [ 0  1  0  0 | v  v  v | v ]
///  pz [ 0  0  1  0 | v  v  v | v ]
///   t [ 0  0  0  1 | 0  0  0 | 0 ]
///  dx [ 0  0  0  0 | v  v  v | v ]
///  dy [ 0  0  0  0 | v  v  v | v ]
///  dz [ 0  0  0  0 | v  v  v | v ]
/// qop [ 0  0  0  0 | 0  0  0 | v ]
///
/// @{
template <bool has_field_gradient>
struct transport_jacobian_substructure_type;

/// No field gradient: the position block stays the identity it started as.
template <>
struct transport_jacobian_substructure_type<false> {
  // clang-format off
  using type = substructure<
      row< one, zero, zero, zero, variable, variable, variable, variable>,
      row<zero,  one, zero, zero, variable, variable, variable, variable>,
      row<zero, zero,  one, zero, variable, variable, variable, variable>,
      row<zero, zero, zero,  one,     zero,     zero,     zero,     zero>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero,     zero,     zero,     zero, variable>
    >::canonical_type;
  // clang-format on
};

/// Field gradient followed: position and direction pick up d/d(position).
template <>
struct transport_jacobian_substructure_type<true> {
  // clang-format off
  using type = substructure<
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<    zero,     zero,     zero,  one,     zero,     zero,     zero,     zero>,
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<variable, variable, variable, zero, variable, variable, variable, variable>,
      row<    zero,     zero,     zero, zero,     zero,     zero,     zero, variable>
    >::canonical_type;
  // clang-format on
};

template <bool has_field_gradient>
using transport_jacobian_substructure =
    typename transport_jacobian_substructure_type<has_field_gradient>::type;
/// @}

static_assert(
    std::is_same_v<transport_jacobian_substructure<false>::multiplication_type<
                       transport_jacobian_substructure<false>>,
                   transport_jacobian_substructure<false>>,
    "The transport Jacobian substructure without field gradient must be closed "
    "under multiplication");

static_assert(
    std::is_same_v<transport_jacobian_substructure<true>::multiplication_type<
                       transport_jacobian_substructure<true>>,
                   transport_jacobian_substructure<true>>,
    "The transport Jacobian substructure with field gradient must be closed "
    "under multiplication");

static_assert(transport_jacobian_substructure<false>::num_variables == 25u);
static_assert(transport_jacobian_substructure<true>::num_variables == 43u);

}  // namespace detray::ksm
