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

  template <typename T>
  using covmat_t = typename T::CovMatrix_t;

  METHOD_TRAIT(reference_surface_t, referenceSurface);

  template <typename T>
  using position_t = decltype(std::declval<T>().position());
  template <typename T>
  using momentum_t = decltype(std::declval<T>().momentum());
  template <typename T>
  using charge_t = decltype(std::declval<T>().charge());
  template <typename T>
  using timet = decltype(std::declval<T>().time());
  template <typename T>
  using covariance_t = decltype(std::declval<T>().covariance());
  template <typename T>
  using parameters_t = decltype(std::declval<T>().parameters());

  // clang-format off
    template <typename P>
      struct ParameterConcept {
        constexpr static bool covmat_exists = exists<covmat_t, const P>;
        static_assert(covmat_exists, "Covariance matrix type not found");
        constexpr static bool reference_surface_exists = has_method<const P, const Surface&, reference_surface_t>;
        static_assert(reference_surface_exists, "referenceSurface method not found");
        constexpr static bool position_exists = identical_to<Vector3D, position_t, const P>;
        static_assert(position_exists, "position method not found");        
        constexpr static bool momentum_exists = identical_to<Vector3D, momentum_t, const P>;
        static_assert(momentum_exists, "momentum method not found");
        constexpr static bool charge_exists = identical_to<double, charge_t, const P>;
        static_assert(charge_exists, "charge method not found");
        constexpr static bool time_exists = identical_to<double, timet, const P>;
        static_assert(time_exists, "time method not found");
        constexpr static bool covariance_exists = identical_to<const std::optional<BoundSymMatrix>&, covariance_t, const P>;
        static_assert(covariance_exists, "covariance method not found");
        constexpr static bool parameters_exists = identical_to<BoundVector, parameters_t, const P>;
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