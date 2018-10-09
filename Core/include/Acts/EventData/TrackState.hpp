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

// @brief enum of readability in parameter acces
enum ParametericType : int { predicted = 0, filtered = 1, smoothed = 2 };

// @brief enum of readability in measurement acces
enum MeasurementType : int { uncalibrated = 0, calibrated = 1 };

/// @brief Parameteric part, non-type dependent.
/// It reduces the number of visitor pattern calls
///
/// This is all the information that concerns the
/// the track parameterisation and the jacobian
/// It is enough to to run the track smoothing
template <typename parameters_t, typename jacobian_t>
struct ParametricState
{
  /// The predicted state
  boost::optional<parameters_t> predicted;
  /// The filtered state
  boost::optional<parameters_t> filtered;
  /// The smoothed state
  boost::optional<parameters_t> smoothed;
  /// The transport jacobian matrix
  boost::optional<jacobian_t> jacobian;
  /// The path length along the track - will help sorting
  double pathLength = 0.;
};

/// @class TrackState
///
/// @brief Templated class to hold the track information
/// on a surface along the trajectory
///
/// @tparam identifier_t Type of the identifier
/// @tparam parameters_t Type of the parameters on the surface
/// @tparam parameters_t Type of the jacobian for the transport
/// @tparam params Type list of the measurement type
///
/// @note the Surface is only stored as a pointer, i.e. it is
/// assumed the surface lives longer than the TrackState
template <typename identifier_t,
          typename parameters_t,
          typename jacobian_t,
          ParID_t... params>
class TrackState
{

public:
  /// @brief Nested measurement part, dependent type.
  /// It reduces the nubmer of vistor pattern calls
  ///
  /// This is the measurement and the calibratedMeasurement
  /// (in case the latter is different)
  struct MeasuredState
  {
    /// The optional measurement
    boost::optional<Measurement<identifier_t, params...>> uncalibrated;
    /// The optional calibrabed measurement
    boost::optional<Measurement<identifier_t, params...>> calibrated;
  };

  /// Constructor from (uncalibrated) measurement (move)
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object
  TrackState(Measurement<identifier_t, params...>&& m)
    : surface(&m.referenceSurface())
  {
    measurement.uncalibrated = std::move(m);
  }

  /// Constructor from (uncalibrated) measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object (moved)
  TrackState(const Measurement<identifier_t, params...>& m)
    : surface(&m.referenceSurface())
  {
    measurement.uncalibrated = std::move(m);
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param p The parameters object (moved)
  TrackState(parameters_t p) : surface(&(p.referenceSurface()))
  {
    parametric.predicted = std::move(p);
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

  /// @brief number of Measured parameters, forwarded
  ///
  /// @return number of measured parameters as a static const expression
  static constexpr unsigned int
  size()
  {
    return Measurement<identifier_t, params...>::size();
  }

  /// The surface of this TrackState
  const Surface* surface = nullptr;

  /// The parametric part
  ParametricState<parameters_t, jacobian_t> parametric;

  /// The measurement part
  MeasuredState measurement;
};

/// @brief track state for measurements
template <typename identifier_t,
          typename parameters_t,
          typename jacobian_t,
          ParID_t... params>
using MeasuredTrackState
    = TrackState<identifier_t, parameters_t, jacobian_t, params...>;

/// @brief track state for parametric description
///
/// @todo: investigate if we can move that to dim=0 description (Eigen allows)
template <typename identifier_t, typename parameters_t, typename jacobian_t>
using ParametricTrackState
    = TrackState<identifier_t, parameters_t, jacobian_t, ParDef::eLOC_0>;

/// @brief general type for any possible Measurement
template <typename identifier_t, typename parameters_t, typename jacobian_t>
using VariantTrackState = typename detail::
    trackstate_type_generator<identifier_t, parameters_t, jacobian_t>::type;
}
