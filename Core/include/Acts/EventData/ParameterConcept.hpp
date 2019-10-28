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

  /// Typedef of CovMatrix_t required
  template <typename T>
  using covmat_t = typename T::CovMatrix_t;

  /// The following lines define functions for surface bound parametrisations
  METHOD_TRAIT(reference_surface_t, referenceSurface);

  template <typename T>
  using position_returntype_t = decltype(std::declval<T>().position());
  template <typename T>
  using momentum_returntype_t = decltype(std::declval<T>().momentum());
  template <typename T>
  using charge_returntype_t = decltype(std::declval<T>().charge());
  template <typename T>
  using time_returntype_t = decltype(std::declval<T>().time());
  template <typename T>
  using covariance_returntype_t = decltype(std::declval<T>().covariance());
  template <typename T>
  using parameters_returntype_t = decltype(std::declval<T>().parameters());

  /// The following lines define functions for free parametrisations
  METHOD_TRAIT(position_t, position);
  METHOD_TRAIT(momentum_t, momentum);
  METHOD_TRAIT(charge_t, charge);
  METHOD_TRAIT(timet, time);
  METHOD_TRAIT(covariance_t, covariance);
  METHOD_TRAIT(parameters_t, parameters);

    template <typename P>
      struct ParameterConcept {
		  
		/// Test if covariance matrix data type is defined
        constexpr static bool covmat_exists = exists<covmat_t, const P>;
        static_assert(covmat_exists, "Covariance matrix type not found");
        
        /// Tests related to surface bound parametrisations
        constexpr static bool reference_surface_exists = has_method<const P, const Surface&, reference_surface_t>;
		constexpr static bool surfacebound_position_exists = identical_to<Vector3D, position_returntype_t, const P>;
		constexpr static bool surfacebound_momentum_exists = identical_to<Vector3D, momentum_returntype_t, const P>;
		constexpr static bool surfacebound_charge_exists = identical_to<double, charge_returntype_t, const P>;
		constexpr static bool surfacebound_time_exists = identical_to<double, time_returntype_t, const P>;
		constexpr static bool surfacebound_covariance_exists = identical_to<const std::optional<BoundSymMatrix>&, covariance_returntype_t, const P>;
		constexpr static bool surfacebound_parameters_exists = identical_to<BoundVector, parameters_returntype_t, const P>;
		
		/// Tie all surface bound results together
		constexpr static bool boundValue = require<reference_surface_exists,
													surfacebound_position_exists,
													surfacebound_momentum_exists,
													surfacebound_charge_exists,
													surfacebound_time_exists,
													surfacebound_covariance_exists,
													surfacebound_parameters_exists>;
		
		/// Tests related to free parametrisations											
		constexpr static bool free_position_exists = has_method<const P, Vector3D, position_t>;
		constexpr static bool free_momentum_exists = has_method<const P, Vector3D, momentum_t>;
		constexpr static bool free_charge_exists = has_method<const P, double, charge_t>;
		constexpr static bool free_time_exists = has_method<const P, double, timet>;
		constexpr static bool free_covariance_exists = has_method<const P, const std::optional<FreeSymMatrix>&, covariance_t>;
		constexpr static bool free_parameters_exists = has_method<const P, FreeVector, parameters_t>;
		
		/// Tie all free results together
		constexpr static bool freeValue = require<free_position_exists,						
													free_momentum_exists,
													free_charge_exists,
													free_time_exists,
													free_covariance_exists,
													free_parameters_exists>;

		/// Assert the individual results for easier error handling
		// TODO: The surface getter should be tested but would be useless in this construction
		// static_assert(reference_surface_exists, "referenceSurface method not found");
		static_assert(surfacebound_position_exists || free_position_exists, "position method not found");        
		static_assert(surfacebound_momentum_exists || free_momentum_exists, "momentum method not found");
		static_assert(surfacebound_charge_exists || free_charge_exists, "charge method not found");
		static_assert(surfacebound_time_exists || free_time_exists, "time method not found");
		static_assert(surfacebound_covariance_exists || free_covariance_exists, "covariance method not found");
		static_assert(surfacebound_parameters_exists || free_parameters_exists, "parameters method not found");

		/// Evaluate that everything required exists
        constexpr static bool value = require<covmat_exists, either<boundValue, freeValue>>;
      };
  }  // namespace Parameter
}  // namespace concept

template <typename parameters_t>
constexpr bool ParameterConcept =
    Acts::concept ::Parameter::ParameterConcept<parameters_t>::value;
}  // namespace Acts