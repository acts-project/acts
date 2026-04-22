// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s).
#include <random>

namespace detray {

template <concepts::algebra algebra_t>
struct scattering_helper {
 public:
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  /// @brief Operator to scatter the direction with scattering angle
  ///
  /// @param dir  input direction
  /// @param angle  scattering angle
  /// @param generator random generator
  /// @returns the new direction from random scattering
  template <typename generator_t>
  DETRAY_HOST inline vector3_type operator()(const vector3_type& dir,
                                             const scalar_type angle,
                                             generator_t& generator) const {
    // Generate theta and phi for random scattering
    const scalar_type r_theta{
        angle == scalar_type{0}
            ? 0.f
            : std::normal_distribution<scalar_type>(0.f, angle)(generator)};

    const scalar_type r_phi{std::uniform_real_distribution<scalar_type>(
        -constant<scalar_type>::pi, constant<scalar_type>::pi)(generator)};

    // xaxis of curvilinear plane
    const vector3_type u =
        unit_vectors<vector3_type>().make_curvilinear_unit_u(dir);

    vector3_type new_dir = axis_rotation<algebra_t>(u, r_theta)(dir);
    return axis_rotation<algebra_t>(dir, r_phi)(new_dir);
  }
};

}  // namespace detray
