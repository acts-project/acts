// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/ratio.hpp"

// System include(s)
#include <tuple>

namespace detray {

/// Compile-time material mixture type. The summation of ratios should be equal
/// to one
template <concepts::scalar scalar_t, typename... material_types>
struct mixture
    : public material<scalar_t, typename ratio_sum<
                                    typename material_types::ratio...>::ratio> {
 public:
  using ratio = typename ratio_sum<typename material_types::ratio...>::ratio;

  static_assert(is_ratio_one_v<ratio>,
                "Sumation of ratios should be equal to 1");

  /// Constructor
  constexpr mixture() {
    // Compute effective relative atomic mass
    // Zeff = Ar0 * ratio0 + Ar1 * ratio1 + ...
    auto sum_Ar = [](material_types... M) constexpr -> decltype(auto) {
      return ((M.Ar() * M.fraction()) + ...);
    };

    this->set_Ar(std::apply(sum_Ar, std::tuple<material_types...>()));

    // Compute effective atomic number
    // Zeff = Z0 * ratio0 + Z1 * ratio1 + ...
    auto sum_Z = [](material_types... M) constexpr -> decltype(auto) {
      return ((M.Z() * M.fraction()) + ...);
    };

    this->set_Z(std::apply(sum_Z, std::tuple<material_types...>()));

    // Get averaged mass density
    auto sum_rho = [](material_types... M) constexpr -> decltype(auto) {
      return ((M.mass_density() * M.fraction()) + ...);
    };

    this->set_mass_density(
        std::apply(sum_rho, std::tuple<material_types...>()));

    // Compute effective radiation length (X0)
    // reference:
    // https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
    // X_avg = W_avg / Sum_i[ W_i / X_i ],
    // where:
    // W_i is mass density of i_th component
    // W_avg is the averaged mass density
    auto sum_rho_over_X0 = [](material_types... M) constexpr -> decltype(auto) {
      return ((M.fraction() / M.X0()) + ...);
    };
    this->set_X0(1.f /
                 std::apply(sum_rho_over_X0, std::tuple<material_types...>()));

    // Compute effective nuclear radiation length
    // Follow the same equation of effective X0
    auto sum_rho_over_L0 = [](material_types... M) constexpr -> decltype(auto) {
      return ((M.fraction() / M.L0()) + ...);
    };

    this->set_L0(1.f /
                 std::apply(sum_rho_over_L0, std::tuple<material_types...>()));

    // Compute molar density
    this->set_molar_density(
        this->mass_to_molar_density(this->Ar(), this->mass_density()));

    // @TODO: Calculate density effect data as well if exist?
  }
};

}  // namespace detray
