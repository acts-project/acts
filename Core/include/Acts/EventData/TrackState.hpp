// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/detail/trackstate_type_generator.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

class Surface;

/// @class TrackState
///
/// @brief Templated class to hold the track information
/// on a surface along the trajectory
///
/// @tparam identifier_t Type of the identifier
/// @tparam parameters_t Type of the parameters on the surface
/// @tparam params Type list of the measurement type
template <typename identifier_t, typename parameters_t, ParID_t... params>
class TrackState
{
public:
  /// Constructor from measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param measurement the object
  TrackState(const Measurement<identifier_t, params...>& measurement)
    : m_surface(&measurement.referenceSurface())
  {
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param parameters object as unitue ptr
  TrackState(std::unique_ptr<const parameters_t> parameters)
    : m_predicted(std::move(parameters))
    , m_surface(&(m_predicted->referenceSurface()))
  {
  }

  /// Constructor from surface
  ///
  /// @tparam parameters_t Type of the predicted parameters
  TrackState(const Surface& surface) : m_surface(&(surface)) {}

  /// @brief return method for the surface
  const Surface&
  referenceSurface() const
  {
    return (*m_surface);
  }

private:
  /// The surface of this TrackState
  const Surface* m_surface = nullptr;

  /// The predicted state if needed
  std::unique_ptr<const parameters_t> m_predicted = nullptr;

  /// The updated state if needed
  std::unique_ptr<const parameters_t> m_updated = nullptr;

  /// The smoothed state if needed
  std::unique_ptr<const parameters_t> m_smoothed = nullptr;

  /// The measurement_t at this TrackState
  std::unique_ptr<const Measurement<identifier_t, params...>> m_measurement
      = nullptr;

};

/// @brief track state for measurements
template <typename identifier_t, typename parameters_t, ParID_t... params>
using MeasuredTrackState = TrackState<identifier_t, parameters_t, params...>;

/// @brief track state for parametric description
///
/// @todo: investigate if we can move that to dim=0 description (Eigen allows)
template <typename identifier_t, typename parameters_t>
using ParametricTrackState
    = TrackState<identifier_t, parameters_t, ParDef::eLOC_0>;

/// @brief general type for any possible Measurement
template <typename identifier_t, typename parameters_t>
using VariantTrackState =
    typename detail::trackstate_type_generator<identifier_t,
                                               parameters_t>::type;
}
