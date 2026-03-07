// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"

namespace Acts::detail {

/// @brief Shorthand for Bound or Free track parameters
template <class parameters_t>
concept isBoundOrFreeTrackParams = FreeTrackParametersConcept<parameters_t> ||
                                   BoundTrackParametersConcept<parameters_t>;

/// @brief Shorthand for GenericBoundTrackParameters
template <class parameters_t>
concept isBoundTrackParams = std::same_as<parameters_t, BoundTrackParameters>;

/// @brief Concept that restricts the type of the
/// accumulation grid cell
template <typename grid_t>
concept TrackParamsGrid = requires {
  typename grid_t::value_type::first_type;
  typename grid_t::value_type::second_type;

  requires isBoundOrFreeTrackParams<
      typename grid_t::value_type::first_type::element_type>;
  requires isBoundOrFreeTrackParams<
      typename grid_t::value_type::second_type::element_type>;

  requires requires(typename grid_t::value_type val) {
    {
      val.first
    } -> std::same_as<
          std::shared_ptr<typename decltype(val.first)::element_type>&>;
    { val.second } -> std::same_as<decltype(val.first)&>;
  };
};

}  // namespace Acts::detail
