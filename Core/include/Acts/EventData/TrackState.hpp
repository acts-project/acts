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
#include "Acts/EventData/detail/trackstate_type_generator.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

class Surface;

// @brief enum of readability in parameter acces
enum class ParametricType : int { predicted = 0, filtered = 1, smoothed = 2 };

// @brief enum of readability in measurement acces
enum class MeasurementType : int { uncalibrated = 0, calibrated = 1 };

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
  boost::optional<parameters_t> predicted{boost::none};
  /// The filtered state
  boost::optional<parameters_t> filtered{boost::none};
  /// The smoothed state
  boost::optional<parameters_t> smoothed{boost::none};
  /// The transport jacobian matrix
  boost::optional<jacobian_t> jacobian{boost::none};
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
/// @tparam jacobian_t Type of the jacobian for the transport
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
    /// The default measured state initializes
    MeasuredState() = default;

    /// The copy constructor
    ///
    /// @param ms The source MeasuredState
    MeasuredState(const MeasuredState& ms)
      : uncalibrated(ms.uncalibrated), calibrated(ms.calibrated)
    {
    }

    /// The move constructor
    ///
    /// @param ms The source MeasuredState
    MeasuredState(MeasuredState&& ms)
      : uncalibrated(std::move(ms.uncalibrated))
      , calibrated(std::move(ms.calibrated))
    {
    }

    /// The copy assignement operator
    ///
    /// @param ms The source MeasuredState
    MeasuredState&
    operator=(const MeasuredState& ms)
    {
      if (this != &ms) {
        uncalibrated = ms.uncalibrated;
        calibrated   = ms.calibrated;
      }
      return (*this);
    }

    /// The Move assignement operator
    ///
    /// @param ms The source MeasuredState
    MeasuredState&
    operator=(MeasuredState&& ms)
    {
      uncalibrated = std::move(ms.uncalibrated);
      calibrated   = std::move(ms.calibrated);
      return (*this);
    }

    /// The optional measurement
    boost::optional<Measurement<identifier_t, params...>> uncalibrated{
        boost::none};
    /// The optional calibrabed measurement
    boost::optional<Measurement<identifier_t, params...>> calibrated{
        boost::none};
  };

  /// Constructor from (uncalibrated) measurement
  ///
  /// @tparam measurement_t Type of the measurement
  /// @param m The measurement object
  TrackState(Measurement<identifier_t, params...> m)
  {
    measurement.uncalibrated = std::move(m);
    assignSurface(measurement.uncalibrated);
  }

  /// Constructor from parameters
  ///
  /// @tparam parameters_t Type of the predicted parameters
  /// @param p The parameters object
  TrackState(parameters_t p)
  {
    parametric.predicted = std::move(p);
    assignSurface(parametric.predicted);
  }

  /// Virtual destructor
  virtual ~TrackState() = default;

  /// Copy constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(const TrackState& rhs)
    : surface(rhs.surface)
    , parametric(rhs.parametric)
    , measurement(rhs.measurement)
  {
  }

  /// Copy move constructor
  ///
  /// @param rhs is the source TrackState
  TrackState(TrackState&& rhs)
    : surface(std::move(rhs.surface))
    , parametric(std::move(rhs.parametric))
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
    parametric  = rhs.parametric;
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
    parametric  = std::move(rhs.parametric);
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
    return Measurement<identifier_t, params...>::size();
  }

  /// The surface of this TrackState
  const Surface* surface = nullptr;

  /// The parametric part
  ParametricState<parameters_t, jacobian_t> parametric;

  /// The measurement part
  MeasuredState measurement;

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
