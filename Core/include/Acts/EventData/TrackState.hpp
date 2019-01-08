// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/none.hpp>
#include <boost/optional.hpp>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
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
template <typename identifier_t, typename parameters_t>
class TrackState
{

public:
  using Identifier = identifier_t;
  using Parameters = parameters_t;
  using Jacobian   = typename Parameters::CovMatrix_t;

  /// Constructor from (uncalibrated) measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object
  TrackState(FittableMeasurement<identifier_t> m)
  {
    this->surface            = MeasurementHelpers::getSurface(m);
    measurement.uncalibrated = std::move(m);
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param p The parameters object
  TrackState(parameters_t p)
  {
    surface             = &p.referenceSurface();
    parameter.predicted = std::move(p);
  }

  /// Virtual destructor
  virtual ~TrackState() = default;

  /// Copy constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(const TrackState& rhs)
    : surface(rhs.surface)
    , parameter(rhs.parameter)
    , measurement(rhs.measurement)
  {
  }

  /// Copy move constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(TrackState&& rhs)
    : surface(std::move(rhs.surface))
    , parameter(std::move(rhs.parameter))
    , measurement(std::move(rhs.measurement))
  {
  }

  /// Assignment operator
  ///
  /// @param rhs is the source TrackState
  TrackState&
  operator=(const TrackState& rhs)
  {
    surface     = rhs.surface;
    parameter   = rhs.parameter;
    measurement = rhs.measurement;
    return (*this);
  }

  /// Assignment move operator
  ///
  /// @param rhs is the source TrackState
  TrackState&
  operator=(TrackState&& rhs)
  {
    surface     = std::move(rhs.surface);
    parameter   = std::move(rhs.parameter);
    measurement = std::move(rhs.measurement);
    return (*this);
  }

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
    return FittableMeasurement<identifier_t>::size();
  }

  /// The surface of this TrackState
  const Surface* surface = nullptr;

  /// The parameter part
  /// This is all the information that concerns the
  /// the track parameterisation and the jacobian
  /// It is enough to to run the track smoothing
  struct
  {
    /// The predicted state
    boost::optional<Parameters> predicted{boost::none};
    /// The filtered state
    boost::optional<Parameters> filtered{boost::none};
    /// The smoothed state
    boost::optional<Parameters> smoothed{boost::none};
    /// The transport jacobian matrix
    boost::optional<Jacobian> jacobian{boost::none};
    /// The path length along the track - will help sorting
    double pathLength = 0.;
  } parameter;

  /// @brief Nested measurement part
  /// This is the uncalibrated and calibrated measurement
  /// (in case the latter is different)
  struct
  {
    /// The optional (uncalibrated) measurement
    boost::optional<FittableMeasurement<identifier_t>> uncalibrated{
        boost::none};
    /// The optional calibrabed measurement
    boost::optional<FittableMeasurement<identifier_t>> calibrated{boost::none};
  } measurement;

private:
  /// Assign the surface from an optional parameter
  ///
  /// @tparam Type of the optional parameter
  template <typename optional_type_t>
  void
  assignSurface(const optional_type_t& optional)
  {
    surface = &(optional.template get().referenceSurface());
  }
};
}
