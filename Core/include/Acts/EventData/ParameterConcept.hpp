// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
class Surface;

namespace concept {
  namespace Parameter {

  METHOD_TRAIT(reference_surface_t, referenceSurface);
  METHOD_TRAIT(position_t, position);
  METHOD_TRAIT(momentum_t, momentum);
  METHOD_TRAIT(charge_t, charge);
  METHOD_TRAIT(timet, time);
  METHOD_TRAIT(covariance_t, covariance);
  METHOD_TRAIT(parameters_t, parameters);

  // clang-format off
    template <typename P>
      struct ParameterConcept {
        constexpr static bool reference_surface_exists = has_method<const P, const Surface&, reference_surface_t>;
        static_assert(reference_surface_exists, "referenceSurface method not found");
        constexpr static bool position_exists = has_method<const P, Vector3D, position_t>;
        static_assert(position_exists, "position method not found");
        constexpr static bool momentum_exists = has_method<const P, Vector3D, momentum_t>;
        static_assert(momentum_exists, "momentum method not found");
        constexpr static bool charge_exists = has_method<const P, double, charge_t>;
        static_assert(charge_exists, "charge method not found");
        constexpr static bool time_exists = has_method<const P, double, timet>;
        static_assert(time_exists, "time method not found");
        constexpr static bool covariance_exists = has_method<const P, const std::optional<BoundSymMatrix>&, covariance_t>;
        static_assert(covariance_exists, "covariance method not found");
        constexpr static bool parameters_exists = has_method<const P, BoundVector, parameters_t>;
        static_assert(parameters_exists, "parameters method not found");
        
        constexpr static bool value = require<reference_surface_exists,
                                              position_exists,
                                              momentum_exists,
                                              charge_exists,
                                              time_exists,
                                              covariance_exists,
                                              parameters_exists>;
      };
  // clang-format on
  }  // namespace Parameter
}  // namespace concept

template <typename parameters_t>
constexpr bool ParameterConcept =
    Acts::concept ::Parameter::ParameterConcept<parameters_t>::value;
}  // namespace Acts
