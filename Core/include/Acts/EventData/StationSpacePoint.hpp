// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/PointerTraits.hpp"

#include <type_traits>

namespace Acts::Experimental {
/// @brief Concept definition of the station space points. They're primarily used in composite detectors,
///        like the Muon stations in side the ATLAS experiment. The stations
///        usually consist of few layers of straw tubes which maybe sandwiched
///        by strip detector layers. The straws are used to measure the passage
///        of the particle in the bending plane, while the strip may supplement
///        the track measurement by providing the measurements along the straw.
template <typename SpacePointType>
concept StationSpacePoint = requires(const SpacePointType sp) {
  ///  @brief Local position of the space point measurement. It's either the position of the wire
  ///         or the position of the fired strip in the station
  { sp.localPosition() } -> std::same_as<const Vector3&>;
  /// @brief Orientation of the sensor, which is either the wire orientation or the strip orientation.
  ///        Distortions along the sensor direction do not alter the track
  ///        residual
  { sp.sensorDirection() } -> std::same_as<const Vector3&>;
  /// @brief Unit vector pointing to the next strip/straw in the plane or
  ///        in case of a combined measurement, the complementary strip
  ///        direction
  { sp.sensorNormal() } -> std::same_as<const Vector3&>;
  /// @brief Normal vector on the strip-plane.
  { sp.planeNormal() } -> std::same_as<const Vector3&>;
  /// @brief Radius of the straw-tube measurement. The returned value is zero for strip measurements
  { sp.driftRadius() } -> std::same_as<double>;
  /// @brief Recorded time of the measurement, if provided by the technology
  { sp.time() } -> std::same_as<double>;
  /// @brief Measurement's covariance. The first two components correspond to
  ///        the local covariances along the non-bending & bending direction
  ///        The third component corresponds to the measurement's time
  ///        covariance.
  { sp.covariance() } -> std::same_as<const std::array<double, 3>&>;

  /// @brief Return whether the space point represents a straw measurement
  { sp.isStraw() } -> std::same_as<bool>;
  /// @brief Return whether the space point provides a direct time constraint
  { sp.hasTime() } -> std::same_as<bool>;
  /// @brief Return whether the space point constrains the bending direction
  ///        of the segment
  { sp.measPrecCoord() } -> std::same_as<bool>;
  /// @brief Return whether the space point constrains the non-bending direction
  ///        of the segment
  { sp.measNonPrecCoord() } -> std::same_as<bool>;
};

/// @brief Define the Space Point pointer concept as an ordinary / smart pointer
///        over space points
template <typename SpacePoint_t>
concept StationSpacePointPtr =
    PointerConcept<SpacePoint_t> &&
    StationSpacePoint<typename RemovePointer<SpacePoint_t>::type>;

/// @brief A station space point container is any std::container over space points
template <typename ContType_t>
concept StationSpacePointContainer =
    requires(ContType_t mCont, const ContType_t cCont) {
      { mCont.begin() } -> std::same_as<typename ContType_t::iterator>;
      { mCont.end() } -> std::same_as<typename ContType_t::iterator>;
      { cCont.begin() } -> std::same_as<typename ContType_t::const_iterator>;
      { cCont.end() } -> std::same_as<typename ContType_t::const_iterator>;
      { cCont.size() } -> std::same_as<typename ContType_t::size_type>;
      { cCont.empty() } -> std::same_as<bool>;
      requires StationSpacePointPtr<typename ContType_t::value_type>;
    };

}  // namespace Acts::Experimental
