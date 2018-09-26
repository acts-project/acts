// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/optional.hpp>
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
///
/// @note the Surface is only stored as a pointer, i.e. it is
/// assumed the surface lives longer than the TrackState
template <typename identifier_t, typename parameters_t, ParID_t... params>
class TrackState
{
public:
  /// Constructor from measurement (move)
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object
  TrackState(Measurement<identifier_t, params...>&& m)
    : surface(&m.referenceSurface()), measurement(std::move(m))
  {
  }

  /// Constructor from measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object
  TrackState(const Measurement<identifier_t, params...>& m)
    : surface(&m.referenceSurface()), measurement(m)
  {
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param p The parameters object as unitue ptr
  TrackState(parameters_t p)
    : surface(&(p.referenceSurface())), predictedState(std::move(p))
  {
  }

  /// Virtual destructor
  virtual ~TrackState() = default;

  /// Copy constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(const TrackState& rhs) = default;

  /// Copy move constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(TrackState&& rhs) = default;

  /// Assignment operator
  ///
  /// @param rhs is the source TrackState
  TrackState&
  operator=(const TrackState& rhs)
      = default;

  /// Assignment move operator
  ///
  /// @param rhs is the source TrackState
  TrackState&
  operator=(TrackState&& rhs)
      = default;

  /// @brief return method for the surface
  const Surface&
  referenceSurface() const
  {
    return (*surface);
  }

  /// The surface of this TrackState
  const Surface* surface = nullptr;

  /// The optional measurement
  boost::optional<Measurement<identifier_t, params...>> measurement;

  /// The optional calibrabed measurement
  boost::optional<Measurement<identifier_t, params...>> calibratedMeasurement;

  /// The predicted state
  boost::optional<parameters_t> predictedState;

  /// The updated state
  boost::optional<parameters_t> updatedState;

  /// The smoothed state
  boost::optional<parameters_t> smoothedState;
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
