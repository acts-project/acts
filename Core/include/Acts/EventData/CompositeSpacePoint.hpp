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
///        of the particle in the bending plane, while the strips may supplement
///        the track measurement by providing coordinates along the straw. The
///        CompositeSpacePoint assumes an orthogonal coordinate system, where
///           x-axis: Is parallel to the straw wires
///           y-axis: Points to the next straw in a layer
///           z-axis: Points outwards from the experiment
template <typename SpacePointType>
concept CompositeSpacePoint = requires(const SpacePointType sp) {
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
  { sp.toNextSensor() } -> std::same_as<const Vector3&>;
  /// @brief Normal vector on the strip-plane spanned by `sensorDirection()` & `toNextSensor()`.
  { sp.planeNormal() } -> std::same_as<const Vector3&>;
  /// @brief Radius of the straw-tube measurement. The returned value is zero for strip measurements
  { sp.driftRadius() } -> std::same_as<double>;
  /// @brief Recorded time of the measurement, if provided by the technology
  { sp.time() } -> std::same_as<double>;
  /// @brief Measurement covariance array. It's composed of the individual strip covariances, making up
  ///        the composite space point and the covariance on time, if provided
  { sp.covariance() } -> std::same_as<const std::array<double, 3>&>;

  /// @brief Return whether the space point represents a straw measurement
  { sp.isStraw() } -> std::same_as<bool>;
  /// @brief Return whether the space point provides a direct time constraint
  { sp.hasTime() } -> std::same_as<bool>;
  /// @brief Returns whether the station space point measures the 0-th coordinate
  ///        and hence constrains the track parameters in the non-bending
  ///        direction
  { sp.measuresLoc0() } -> std::same_as<bool>;
  /// @brief Returns whether the station space point measures the 1-st coordinate
  ///        and hence constrains the track parameters in the bending direction
  { sp.measuresLoc1() } -> std::same_as<bool>;
  /// @brief Returns the spatial dimension of the space point
  { sp.dimension() } -> std::same_as<unsigned>;
};

/// @brief Define the space point pointer concept as an ordinary / smart pointer
///        over space points
template <typename SpacePoint_t>
concept CompositeSpacePointPtr =
    PointerConcept<SpacePoint_t> &&
    CompositeSpacePoint<typename RemovePointer<SpacePoint_t>::type>;

/// @brief A station space point container is any std::container over space points
template <typename ContType_t>
concept CompositeSpacePointContainer =
    std::move_constructible<ContType_t> &&
    requires(ContType_t mCont, const ContType_t cCont) {
      { mCont.begin() } -> std::same_as<typename ContType_t::iterator>;
      { mCont.end() } -> std::same_as<typename ContType_t::iterator>;
      { cCont.begin() } -> std::same_as<typename ContType_t::const_iterator>;
      { cCont.end() } -> std::same_as<typename ContType_t::const_iterator>;
      { cCont.size() } -> std::same_as<typename ContType_t::size_type>;
      { cCont.empty() } -> std::same_as<bool>;
      requires CompositeSpacePointPtr<typename ContType_t::value_type>;
    };

}  // namespace Acts::Experimental
