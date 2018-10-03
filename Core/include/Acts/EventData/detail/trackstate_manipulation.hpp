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
    /// @pram psType Type of the TrackState
    explicit ParametersGetter(StateType psType) : sType(psType) {}

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
      case predicted:
        return edm.predictedState;
      case updated:
        return edm.updatedState;
      default:
        return edm.smoothedState;
      }
    }
    /// The state type for the retrieving
    StateType sType = StateType::predicted;
  };

  /// @brief Visitor pattern to extract the surface
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  template <typename parameters_t>
  struct ParametersSetter : public boost::static_visitor<void>
  {
  public:
    /// Explicit constructor
    ///
    /// @pram psType Type of the TrackState
    explicit ParametersSetter(parameters_t pars, StateType psType)
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
      case predicted:
        edm.predictedState = std::move(sParameters);
        break;
      case updated:
        edm.updatedState = std::move(sParameters);
        break;
      default:
        edm.smoothedState = std::move(sParameters);
      }
    }
    /// The parameters that will be moved into the track state
    parameters_t sParameters;

    /// The type of state that will be updated
    StateType sType = StateType::predicted;
  };

  /// @brief get method to be used with the visitor pattern
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename parameters_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  const boost::optional<parameters_t>&
  getParamaters(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm,
                StateType                                           sType)
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
                StateType                                     sType)
  {
    ParametersSetter<parameters_t> ps(std::move(pars), sType);
    return boost::apply_visitor(ps, edm);
  }

}  // namespace detail
}  // namespace Acts