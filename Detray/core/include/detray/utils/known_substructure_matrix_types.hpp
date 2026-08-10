// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/utils/known_substructure_matrix.hpp"

// System include(s)
#include <type_traits>

/// @file
/// Substructures of the matrices detray actually computes, as opposed to the
/// generic building blocks in known_substructure_matrix.hpp.

namespace detray::ksm {

/// @brief Substructure of the full Jacobian between two bound track states.
///
/// The Jacobian that @c update_full_jacobian_*_impl produces is not dense.
/// Two rows are special, for quite different reasons.
///
/// Bound time decouples unconditionally. The free-to-bound row for time is the
/// unit vector on free time, the path-to-free derivative for free time is hard
/// zero, and the transport Jacobian's time row is a unit vector too, so the
/// product leaves row 5 as e5 and column 5 as zero above the diagonal. Nothing
/// in the configuration changes that.
///
/// Bound q/p decouples only when the step carries no volume material. Then
/// @c dqopds is exactly @c 0.f and @c d2qopdsdqop is too -- so @c dqopqop is
/// exactly one -- and row 4 collapses to e4. With volume material the q/p row
/// picks up the energy-loss couplings and is fully free. Note that this is
/// governed by the presence of volume material alone: @c use_eloss_gradient
/// only forces @c dqopqop to one, which the no-material case already gives.
///
/// The q/p *column* stays free either way: momentum influences the track
/// through curvature, which is a magnetic-field effect and independent of
/// material.
///
///     has_volume_material = true          has_volume_material = false
///        l0  l1 phi  th q/p   t              l0  l1 phi  th q/p   t
///  l0 [[  v,  v,  v,  v,  v,  0],      l0 [[  v,  v,  v,  v,  v,  0],
///  l1  [  v,  v,  v,  v,  v,  0],      l1  [  v,  v,  v,  v,  v,  0],
/// phi  [  v,  v,  v,  v,  v,  0],     phi  [  v,  v,  v,  v,  v,  0],
///  th  [  v,  v,  v,  v,  v,  0],      th  [  v,  v,  v,  v,  v,  0],
/// q/p  [  v,  v,  v,  v,  v,  0],     q/p  [  0,  0,  0,  0,  1,  0],
///   t  [  0,  0,  0,  0,  0,  1]]       t  [  0,  0,  0,  0,  0,  1]]
///
///          25 of 36 free                        20 of 36 free
///
/// @c true is the weaker claim and holds in both cases, so it is the safe
/// choice where the presence of volume material is not known statically.
/// @{
template <bool has_volume_material>
struct full_jacobian_substructure_type;

/// Volume material present: the q/p row carries the energy-loss couplings.
template <>
struct full_jacobian_substructure_type<true> {
  using type =
      substructure<row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<zero, zero, zero, zero, zero, one>>::canonical_type;
};

/// No volume material: q/p is unchanged along the step, so its row is e4.
template <>
struct full_jacobian_substructure_type<false> {
  using type =
      substructure<row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<variable, variable, variable, variable, variable, zero>,
                   row<zero, zero, zero, zero, one, zero>,
                   row<zero, zero, zero, zero, zero, one>>::canonical_type;
};

template <bool has_volume_material>
using full_jacobian_substructure =
    typename full_jacobian_substructure_type<has_volume_material>::type;
/// @}

/// Both variants are closed under multiplication: the product of two full
/// Jacobians has the same substructure again. That is what lets an accumulated
/// Jacobian keep its structure, and its storage, over arbitrarily many steps
/// rather than degrading to dense after the first product.
/// @{
static_assert(
    std::is_same_v<full_jacobian_substructure<true>::multiplication_type<
                       full_jacobian_substructure<true>>,
                   full_jacobian_substructure<true>>,
    "the full Jacobian substructure with volume material is not closed under "
    "multiplication: an accumulated Jacobian would not keep it");

static_assert(
    std::is_same_v<full_jacobian_substructure<false>::multiplication_type<
                       full_jacobian_substructure<false>>,
                   full_jacobian_substructure<false>>,
    "the full Jacobian substructure without volume material is not closed "
    "under multiplication: an accumulated Jacobian would not keep it");
/// @}

/// @brief Substructure of the RK transport Jacobian.
///
/// The 8x8 free-to-free Jacobian the Runge-Kutta stepper accumulates. Both the
/// step matrix D and the running Jacobian carry it, and it is the structure
/// that @c gen_transport_jacobian_types.py generates a bespoke storage type
/// for. Free time is untouched by the step, so its row is a unit vector; the
/// position rows couple to position only when the stepper follows the spatial
/// gradient of the magnetic field.
///
///     has_field_gradient = false            has_field_gradient = true
///     [ 1  0  0  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  1  0  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  0  1  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  0  0  1 | 0  0  0 | 0 ]          [ 0  0  0  1 | 0  0  0 | 0 ]
///     [ 0  0  0  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  0  0  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  0  0  0 | v  v  v | v ]          [ v  v  v  0 | v  v  v | v ]
///     [ 0  0  0  0 | 0  0  0 | v ]          [ 0  0  0  0 | 0  0  0 | v ]
///
///           25 of 64 free                         43 of 64 free
/// @{
template <bool has_field_gradient>
struct transport_jacobian_substructure_type;

/// No field gradient: the position block stays the identity it started as.
template <>
struct transport_jacobian_substructure_type<false> {
  using type = substructure<
      row<one, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, one, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, one, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, one, zero, zero, zero, zero>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero, variable, variable, variable, variable>,
      row<zero, zero, zero, zero, zero, zero, zero, variable>>::canonical_type;
};

/// Field gradient followed: position and direction pick up d/d(position).
template <>
struct transport_jacobian_substructure_type<true> {
  using type = substructure<
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<zero, zero, zero, one, zero, zero, zero, zero>,
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<variable, variable, variable, zero, variable, variable, variable,
          variable>,
      row<zero, zero, zero, zero, zero, zero, zero, variable>>::canonical_type;
};

template <bool has_field_gradient>
using transport_jacobian_substructure =
    typename transport_jacobian_substructure_type<has_field_gradient>::type;
/// @}

/// Closed under multiplication, so the Jacobian the stepper accumulates step
/// by step keeps its structure rather than degrading to dense.
/// @{
static_assert(
    std::is_same_v<transport_jacobian_substructure<false>::multiplication_type<
                       transport_jacobian_substructure<false>>,
                   transport_jacobian_substructure<false>>,
    "the transport Jacobian substructure without field gradient is not closed "
    "under multiplication: the accumulated Jacobian would not keep it");

static_assert(
    std::is_same_v<transport_jacobian_substructure<true>::multiplication_type<
                       transport_jacobian_substructure<true>>,
                   transport_jacobian_substructure<true>>,
    "the transport Jacobian substructure with field gradient is not closed "
    "under multiplication: the accumulated Jacobian would not keep it");

static_assert(transport_jacobian_substructure<false>::num_variables == 25u);
static_assert(transport_jacobian_substructure<true>::num_variables == 43u);
/// @}

}  // namespace detray::ksm
