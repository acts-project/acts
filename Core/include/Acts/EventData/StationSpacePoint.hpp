// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/PointerConcept.hpp"

#include <type_traits>

namespace Acts {
/** @brief Concept definition of the station space points. They're primarly used in composite detectors,
 *         like the Muon stations in side the ATLAS experiment. The stations
 * usually consist of few layers of straw tubes which maybe sandwiched by other
 * strip detector layers. The straws are used to measure the passage of the
 * particle in the bending plane, while the strip may supplement the track
 * measurement by providing measurements along the straw.
 *
 */
template <typename SpacePointType>
concept StationSpacePoint = requires(const SpacePointType sp) {
  /** @brief Local position of the space point measurement. It'either
   *         the position of the wire or the position of the fired strip in the
   * station */
  { sp.localPosition() } -> std::same_as<const Acts::Vector3&>;
  /** @brief Orientation of the sensor, which is either the wire orientation or
   *         the strip orientation. Travelling along the direction does not
   * alter the residual */
  { sp.sensorDirection() } -> std::same_as<const Acts::Vector3&>;
  /** @brief Normal vector on the strip-plane. */
  { sp.stripPlaneNormal() } -> std::same_as<const Acts::Vector3&>;
  /** @brief Radius of the  tube drift-measurement. The returned value is zero for strip measurements */
  { sp.driftRadius() } -> std::same_as<double>;
  /** @brief Time when the measurement was taken */
  { sp.time() } -> std::same_as<double>;
  /** @brief Return the space point covariance. The upper left 2x2 block
   *         describes the spatial meaurement uncertainty. The remaining block
   *         the correlation between the time <-> spatial measurement or the
   * pure time resolution */
  { sp.covariance() } -> std::same_as<const ActsSquareMatrix<3>&>;
};

/** @brief Define the Space Point pointer concept as an ordinary / smart pointer
 *         over space points */
template <typename SpacePoint_t>
concept StationSpacePointPtr =
    PointerConcept<SpacePoint_t> &&
    StationSpacePoint<typename std::remove_pointer<SpacePoint_t>::type>;

/** @brief A station space point container is any std::container over space points */
template <typename ContType_t>
concept StationSpacePointContainer =
    requires(ContType_t cont, const ContType_t const_cont) {
      { cont.begin() } -> std::same_as<typename ContType_t::iterator>;
      { cont.end() } -> std::same_as<typename ContType_t::iterator>;
      {
        const_cont.begin()
      } -> std::same_as<typename ContType_t::const_iterator>;
      { const_cont.end() } -> std::same_as<typename ContType_t::const_iterator>;
      { const_cont.size() } -> std::same_as<typename ContType_t::size_type>;
      { const_cont.empty() } -> std::same_as<bool>;
      requires StationSpacePointPtr<typename ContType_t::value_type>;
    };

}  // namespace Acts