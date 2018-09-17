// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

class Surface;

/// @class TrackState
///
/// @brief Templated class to hold the track information
/// on a surface along the trajectory
///
/// @tparam parameters_t Type of the parameters on the surface
/// @tparam states_t Type of the state object
template <typename parameters_t, typename measurement_t>
class TrackState
{
public:
  /// The surface of this TrackState
  const Surface* surface = nullptr;

  /// The predicted state if needed
  std::unique_ptr<const parameters_t> predicted = nullptr;

  /// The updated state if needed
  std::unique_ptr<const parameters_t> updated = nullptr;

  /// The smoothed state if needed
  std::unique_ptr<const parameters_t> smoothed = nullptr;

  /// The measurement_t at this TrackState
  std::unique_ptr<const measurement_t> measurement = nullptr;

  /// Material Properties associated to this TrackState
  MaterialProperties material{};

  /// Constructor from measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param measurement the object
  TrackState(const measurement_t& measurement)
    : surface(&measurement.referenceSurface())
  {
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param parameters object as unitue ptr
  TrackState(std::unique_ptr<const parameters_t> parameters)
    : surface(&(parameters->referenceSurface()))
    , predicted(std::move(parameters))
  {
  }
};
}
