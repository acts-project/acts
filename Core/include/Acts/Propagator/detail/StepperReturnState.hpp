// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace detail {

	/// @brief This is a helper struct to deduce the dimensions of the full Jacobian. It decides based on the start and end parameters which one it will be. This leads to four different cases ...
	///
	/// @tparam S The boolean expression whether the start parameters are in local representation
	/// @tparam E The boolean expression whether the end parameters are in local representation
  template<bool S, bool E>
  struct JacobianHelper;

  /// @brief Case: 
  /// - start parameters are in local representation
  /// - end parameters are in local representation
  template<>
  struct JacobianHelper<true, true>
  {
	  using type = BoundMatrix;
  };
  
  /// @brief Case: 
  /// - start parameters are in local representation
  /// - end parameters are in global representation  
  template<>
  struct JacobianHelper<true, false>
  {
	  using type = BoundToFreeMatrix;
  };
 
   /// @brief Case: 
  /// - start parameters are in global representation
  /// - end parameters are in local representation 
template<>
  struct JacobianHelper<false, true>  
  {
	  using type = FreeToBoundMatrix;
  };
  
  /// @brief Case: 
  /// - start parameters are in global representation
  /// - end parameters are in global representation
template<>
  struct JacobianHelper<false, false>  
  {
	  using type = FreeMatrix;
  };
  
 // TODO: There should be a deduction of the charge policy
  // This struct is a meta-function deduces the return state if it ends on a surface...
  template <bool local_start, typename end_parameters_t, typename Surface>
  struct s {
	  // The dimensions of the jacobian
	  using JacobianType = typename JacobianHelper<local_start, true>::type;
	// The returning state
	using type = std::tuple<BoundParameters, const JacobianType, double>;
  };

  // ...and for the case that it does not.
  template <bool local_start, typename end_parameters_t>
  struct s<local_start, end_parameters_t, int> {
	  	  // The dimensions of the jacobian
	using JacobianType = typename JacobianHelper<local_start, end_parameters_t::is_local_representation>::type;
		// The returning state
    using type = std::tuple<end_parameters_t, const JacobianType, double>;
  };
  
  /// Return parameter state deduction dependin on the propagation mode
  template <bool local_start, typename end_parameters_t, typename surface_t = int> // TODO: the start parameters are usually unknown or guessed - if jacToGlobal exists then it would be more reliable
  using return_state_type = typename s<local_start, end_parameters_t, surface_t>::type;
  
}  // namespace detail
}  // namespace Acts
