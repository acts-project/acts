// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// boost include(s)
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @cond detail
namespace detail {

  /// @brief Visitor pattern to extract the optional boost parameter
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  template <typename parameters_t>
  struct ParametersGetter
      : public boost::static_visitor<const boost::optional<parameters_t>&>
  {
  public:
    /// Explicit constructor
    ///
    /// @param psType Type of the TrackState
    explicit ParametersGetter(ParametricType psType) : sType(psType) {}

    /// @brief Call operator for parameters extraction using the boost visitor
    /// pattern
    ///
    /// @tparam track_state_t Type of the measurement object (templated)
    /// @param edm The edm object for which the parameters will be extracted
    template <typename track_state_t>
    const boost::optional<parameters_t>&
    operator()(const track_state_t& edm) const
    {
      switch (sType) {
      case ParametricType::predicted:
        return edm.parametric.predicted;
      case ParametricType::filtered:
        return edm.parametric.filtered;
      default:
        return edm.parametric.smoothed;
      }
    }
    /// The state type for the retrieving
    ParametricType sType = ParametricType::predicted;
  };

  /// @brief Visitor pattern to extract the surface
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  template <typename parameters_t>
  struct ParametersSetter : public boost::static_visitor<void>
  {
  public:
    /// Default constructor is deleted
    ParametersSetter() = delete;

    /// Explicit constructor
    ///
    /// @pram psType Type of the TrackState
    explicit ParametersSetter(parameters_t pars, ParametricType psType)
      : sParameters(std::move(pars)), sType(psType)
    {
    }

    /// @brief Call operator for parameters seting using the boost visitor
    /// pattern
    ///
    /// @tparam track_state_t Type of the measurement object (templated)
    /// @param edm The edm object for which the parameters will be extracted
    template <typename track_state_t>
    void
    operator()(track_state_t& edm)
    {
      switch (sType) {
      case ParametricType::predicted:
        edm.parametric.predicted = std::move(sParameters);
        break;
      case ParametricType::filtered:
        edm.parametric.filtered = std::move(sParameters);
        break;
      default:
        edm.parametric.smoothed = std::move(sParameters);
      }
    }
    /// The parameters that will be moved into the track state
    parameters_t sParameters;

    /// The type of state that will be filtered
    ParametricType sType = ParametricType::predicted;
  };

  /// @brief Visitor pattern to extract the optional boost parameter
  ///
  /// @tparam parametric_state_t the bound parametric state
  template <typename parametric_state_t>
  struct ParametricStateGetter
      : public boost::static_visitor<parametric_state_t&>
  {
  public:
    /// Explicit constructor
    explicit ParametricStateGetter() = default;

    /// @brief Call operator for extracting the parameteric_state_t
    ///
    /// @tparam track_state_t Type of the measurement object (templated)
    /// @tparam parametric_state_t Type of the parameteric state
    ///
    /// @param edm The edm object for which the parameters will be extracted
    template <typename track_state_t>
    parametric_state_t&
    operator()(track_state_t& edm) const
    {
      return edm.parametric;
    }
  };

  /// @brief get parametric state with a visitor pattern
  ///
  /// @tparam parametric_state_t Type of the parameteric state
  ///
  /// @param edm The edm object for which the parameters will be extracted
  template <typename parametric_state_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  parametric_state_t&
  getParametricState(boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm)
  {
    static const ParametricStateGetter<parametric_state_t> psg
        = ParametricStateGetter<parametric_state_t>();
    return boost::apply_visitor(psg, edm);
  }

  /// @brief get method to be used with the visitor pattern
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename parameters_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  boost::optional<parameters_t>
  getParamaters(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm,
                ParametricType                                      sType)
  {
    ParametersGetter<parameters_t> pg(sType);
    return boost::apply_visitor(pg, edm);
  }

  /// @brief set method to be used with the visitor pattern
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename parameters_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  void
  setParameters(boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm,
                parameters_t                                  pars,
                ParametricType                                sType)
  {
    ParametersSetter<parameters_t> ps(std::move(pars), sType);
    return boost::apply_visitor(ps, edm);
  }

  /// @brief Visitor pattern to extract the optional boost parameter
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  template <typename measurement_t>
  struct MeasurementGetter
      : public boost::static_visitor<const boost::optional<measurement_t>&>
  {
  public:
    /// Explicit constructor
    ///
    /// @param pmType Type of the Measurement
    explicit MeasurementGetter(MeasurementType pmType) : mType(pmType) {}

    /// @brief Call operator for measurement extraction using the boost visitor
    /// pattern
    ///
    /// @tparam track_state_t Type of the measurement object (templated)
    /// @param edm The edm object for which the measurement will be extracted
    template <typename track_state_t>
    const boost::optional<measurement_t>&
    operator()(const track_state_t& edm) const
    {
      switch (mType) {
      case MeasurementType::uncalibrated:
        return edm.measurement.uncalibrated;
      default:
        return edm.measurement.calibrated;
      }
    }
    /// The state type for the retrieving
    MeasurementType mType = MeasurementType::uncalibrated;
  };

  /// @brief get method to be used with the visitor pattern
  ///
  /// @tparam measurement_t the measurement type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename measurement_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  boost::optional<measurement_t>
  getMeasurement(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm,
                 ParametricType                                      mType)
  {
    MeasurementGetter<measurement_t> mg(mType);
    return boost::apply_visitor(mg, edm);
  }

  /// @ brief visitor pattern to extract the path length
  ///
  struct PathLengthGetter : public boost::static_visitor<double>
  {
  public:
    /// @brief call operator for surface extraction using the boost visitor
    /// pattern
    /// @tparam measurement_t Type of the measurement (templated)
    /// @param m The measurement for which the surface will be extracted
    template <typename edm_object_t>
    double
    operator()(const edm_object_t& edm) const
    {
      return edm.parametric.pathLength;
    }
  };

  template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
  double
  getPathLength(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm)
  {
    static const PathLengthGetter plg = PathLengthGetter();
    return boost::apply_visitor(plg, edm);
  }

}  // namespace detail
}  // namespace Acts