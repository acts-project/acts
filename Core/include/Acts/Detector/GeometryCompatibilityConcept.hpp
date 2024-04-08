// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts::Concepts {

// Types to check compatibility of
template <typename propagator_state_t, typename navigator_t>
using ReturnTypeCurrentVolume =
    decltype(std::declval<navigator_t>().currentVolume(
        std::declval<propagator_state_t>().navigation));

/// @brief Concept ensuring compatibility TrackingGeometry
/// and Detector navigation interfaces with the client code
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam navigator_t Type of the navigator object
template <typename propagator_state_t, typename navigator_t>
struct NavigationCompatibilityConceptImpl {
  /// @brief Ensure that the currentVolume method
  /// returns one of the known volume types
  constexpr static bool isCurrentVolumePtr =
      (Acts::Concepts::identical_to<const TrackingVolume*,
                                    ReturnTypeCurrentVolume, propagator_state_t,
                                    navigator_t> ||
       Acts::Concepts::identical_to<const Acts::Experimental::DetectorVolume*,
                                    ReturnTypeCurrentVolume, propagator_state_t,
                                    navigator_t>);

  static_assert(isCurrentVolumePtr,
                "Return type is not a known volume pointer type");

  constexpr static bool value = Acts::Concepts::require<isCurrentVolumePtr>;
};

template <typename propagator_state_t, typename navigator_t>
constexpr bool NavigationCompatibilityConcept =
    NavigationCompatibilityConceptImpl<propagator_state_t, navigator_t>::value;

}  // namespace Acts::Concepts
